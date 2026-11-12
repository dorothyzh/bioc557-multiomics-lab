/* Analysis worker: keeps the dataset, runs the pipeline on request, returns plain JSON. */
importScripts('engine.js?v=6');
const E = self.Engine;
let DS = null, GS = null, LAST = null;

const send = (type, payload) => postMessage(Object.assign({ type }, payload));
const progress = (id, msg) => send('progress', { id, msg });

self.onmessage = (ev) => {
  const { type, id } = ev.data;
  try {
    if (type === 'load') { DS = ev.data.dataset; GS = ev.data.genesets || null; LAST = null; send('loaded', { id }); return; }
    if (type === 'genesets') { GS = ev.data.genesets; send('done', { id, result: true }); return; }
    if (!DS) throw new Error('No dataset loaded.');
    if (type === 'run') return send('done', { id, result: run(ev.data.params, id) });
    if (type === 'kexp') return send('done', { id, result: kExperiment(ev.data.Ks, id) });
    if (type === 'enrich') return send('done', { id, result: enrich(ev.data.q) });
    if (type === 'network') return send('done', { id, result: network(ev.data.q) });
    if (type === 'predict') return send('done', { id, result: predict(ev.data.q, id) });
    if (type === 'stability') return send('done', { id, result: stability(ev.data.q, id) });
    if (type === 'concordance') return send('done', { id, result: concordance(ev.data.q) });
    throw new Error('Unknown request ' + type);
  } catch (err) { send('error', { id, message: err.message || String(err), stack: err.stack }); }
};

// ---------- sample bookkeeping ----------
function contrastCodes(samples) { const lv = DS.contrast.levels; return samples.map(s => (s[DS.contrast.key] === lv[1] ? 1 : 0)); }

/* Build the per-view sample alignment for the chosen challenge.
   Returns {samples (array of sample objects in analysis order), rows: per view array of row indices into view.values} */
function alignment(params) {
  const rand = E.rng(params.seed || 557);
  let samples = DS.samples.slice();
  if (params.challenge === 'subsample') {
    const frac = params.subsampleFrac || 0.5;
    if (DS.pairing) {
      const ids = [...new Set(samples.map(s => s[DS.pairing.key]))]; const keep = new Set(E.shuffleIdx(ids.length, rand).slice(0, Math.max(4, Math.round(ids.length * frac))).map(i => ids[i]));
      samples = samples.filter(s => keep.has(s[DS.pairing.key]));
    } else {
      const lv = DS.contrast.levels; const out = [];
      lv.forEach(l => { const grp = samples.filter(s => s[DS.contrast.key] === l); E.shuffleIdx(grp.length, rand).slice(0, Math.max(3, Math.round(grp.length * frac))).forEach(i => out.push(grp[i])); });
      const set = new Set(out.map(s => s.id)); samples = samples.filter(s => set.has(s.id));
    }
  }
  // each view stores its own sample_ids; match by ID (never by position)
  const rows = DS.views.map((v, m) => {
    let ids = v.sample_ids.slice(); let order = ids.map((_, i) => i);
    if (params.challenge === 'shuffle') order = E.shuffleIdx(ids.length, E.rng(99 + m)); // rows move together with their IDs
    const pos = new Map(order.map(i => [ids[i], i]));
    let r = samples.map(s => pos.get(s.id));
    if (params.challenge === 'mispair' && m === DS.views.length - 1) { // values move, labels stay: a derangement
      const p = E.shuffleIdx(r.length, E.rng(2024)); for (let i = 0; i < p.length; i++) if (p[i] === i) { const j = (i + 1) % p.length; [p[i], p[j]] = [p[j], p[i]]; }
      r = p.map(i => r[i]);
    }
    return r;
  });
  return { samples, rows };
}

function maskedView(v, frac, seed) { // hide a random fraction of observed values (teaching perturbation)
  const rand = E.rng(seed);
  return Object.assign({}, v, { values: v.values.map(row => row.map(x => (x != null && rand() < frac ? null : x))) });
}

// ---------- main pipeline ----------
function run(params, id) {
  const t0 = Date.now(), T = {}; let tl = t0; const lap = k => { const t = Date.now(); T[k] = t - tl; tl = t; };
  const { samples, rows } = alignment(params);
  const n = samples.length, g = contrastCodes(samples);
  const pairs = DS.pairing && params.paired !== false ? samples.map(s => s[DS.pairing.key]) : null;
  progress(id, 'Preparing each layer');
  const srcViews = params.challenge === 'mask' ? DS.views.map((v, m) => maskedView(v, params.maskFrac || 0.1, 31 + m)) : DS.views;
  const batchCodes = params.batchAdjust && DS.batch ? samples.map(s => s[DS.batch]) : null;
  const prep = srcViews.map((v, m) => E.prepareView(v, rows[m], Object.assign({}, params, { batchCodes })));
  const views = prep.map((p, m) => {
    const v = DS.views[m];
    // QC
    const qc = {
      dInput: p.dInput, dPass: p.keep.length, dModel: p.k, missingTotal: countMiss(p.X) / p.X.length, nImputedForPCA: countMiss(p.M),
      sampleMissing: Array.from({ length: n }, (_, i) => { let c = 0; for (let j = 0; j < p.d; j++) if (p.X[i * p.d + j] !== p.X[i * p.d + j]) c++; return c / p.d; }),
      sampleBox: p.sampleBox,
      featureMean: Array.from(p.featureMean), featureMissing: Array.from({ length: p.d }, (_, j) => { let c = 0; for (let i = 0; i < n; i++) if (p.X[i * p.d + j] !== p.X[i * p.d + j]) c++; return c / n; })
    };
    // sample correlation on model features (pairwise complete, Pearson)
    const corr = sampleCorrelation(p.M, n, p.k);
    // PCA on imputed model matrix
    const Mi = E.imputeColumns(p.M, n, p.k, params.impute || 'median');
    const pc = E.pca(Mi, n, p.k, 5);
    // differential test on all retained features
    const covDesign = !pairs && params.adjust && params.adjust.length ? designFor(samples, params.adjust) : null;
    const de = E.featureTests(p.X, n, p.d, p.keep, g, covDesign, { minPerGroup: 3, pairs });
    const deRows = de.results.map(r => ({ j: r.j, f: v.features[r.j], label: v.labels[r.j], sym: v.symbols[r.j], diff: r.diff, t: r.t, p: r.p, q: r.q, welchP: r.welchP, welchQ: r.welchQ, n0: r.n0, n1: r.n1, nPairs: r.nPairs, mean0: r.mean0, mean1: r.mean1 }));
    return { key: v.key, name: v.name, short: v.short, unit: v.unit, qc, corr, pca: { scores: pc.scores.slice(0, 5), varExp: pc.varExp.slice(0, 5), loadings: pc.loadings }, de: deRows, d0: de.d0, s02: de.s02,
      modelFeatures: p.modelCols.map(j => j), modelLabels: p.modelCols.map(j => v.labels[j]), modelSymbols: p.modelCols.map(j => v.symbols[j]) };
  });
  lap('layers');
  progress(id, 'Joint PCA');
  // joint PCA: each block scaled by 1/sqrt(k) after imputation
  const blocks = prep.map(p => { const Y = E.imputeColumns(p.M, n, p.k, 'median'); const w = 1 / Math.sqrt(p.k); for (let i = 0; i < Y.length; i++) Y[i] *= w; return { Y, k: p.k }; });
  const ktot = blocks.reduce((a, b) => a + b.k, 0), C = new Float64Array(n * ktot); let off = 0;
  blocks.forEach(b => { for (let i = 0; i < n; i++) for (let j = 0; j < b.k; j++) C[i * ktot + off + j] = b.Y[i * b.k + j]; off += b.k; });
  const jp = E.pca(C, n, ktot, 5);
  lap('jointPCA');
  progress(id, 'Fitting the multi-view factor model');
  const fit = E.gfa(prep.map(p => ({ M: p.M, k: p.k })), n, params.K || 5, { seed: params.seed || 557, onIter: (it) => progress(id, `Fitting the multi-view factor model · iteration ${it}`) });
  lap('gfa');
  const weights = fit.W.map((W, m) => W.map((w, j) => w)); // [view][feature][factor]
  // factor-covariate associations
  const covs = DS.covariates.filter(c => c.type !== 'id');
  const assoc = covs.map(c => {
    const use = samples.map(s => !(c.tumour_only && s[c.key] === 'Normal tissue') && s[c.key] != null && s[c.key] !== '');
    return { key: c.key, label: c.label, type: c.type, perFactor: Array.from({ length: fit.K }, (_, k) => {
      const sc = [], codes = []; samples.forEach((s, i) => { if (use[i]) { sc.push(fit.Z[i][k]); codes.push(c.type === 'numeric' ? +s[c.key] : s[c.key]); } });
      const a = E.associate(sc, { numeric: c.type === 'numeric', codes }); return { stat: a.stat, p: a.p, kind: a.kind, n: sc.length };
    }) };
  });
  // view agreement: Spearman correlation between within-view sample distance matrices (Mantel statistic)
  const dist = blocks.map(b => sampleDistances(b.Y, n, b.k));
  const agree = []; for (let a = 0; a < dist.length; a++) for (let b = a + 1; b < dist.length; b++) agree.push({ a, b, r: E.spearman(dist[a], dist[b]).r });
  // variance partition
  lap('assoc+agree');
  progress(id, 'Partitioning variance');
  const vp = prep.map((p, m) => variancePartition(p, samples, g, pairs, params));
  const res = {
    params, n, samples: samples.map(s => s.id), groups: samples.map(s => s[DS.contrast.key]), codes: g, paired: !!pairs,
    views, joint: { scores: jp.scores.slice(0, 5), varExp: jp.varExp.slice(0, 5) },
    gfa: { Z: fit.Z, W: weights, r2: fit.r2, K: fit.K, iterations: fit.iterations, converged: fit.converged },
    assoc, agree, vp, ms: Date.now() - t0, timing: (lap('vp'), T)
  };
  if (!params.challenge || params.challenge === 'none') LAST = { params, prep, samples, g, fit, pairs };
  return res;
}
function countMiss(X) { let c = 0; for (let i = 0; i < X.length; i++) if (X[i] !== X[i]) c++; return c; }
function designFor(samples, keys) {
  const cols = [];
  keys.forEach(k => {
    const c = DS.covariates.find(x => x.key === k); if (!c) return;
    if (c.type === 'numeric') { const v = samples.map(s => +s[k]); const m = E.mean(v); cols.push(v.map(x => (x === x ? x - m : 0))); }
    else { const lv = [...new Set(samples.map(s => s[k]))].sort(); lv.slice(1).forEach(l => cols.push(samples.map(s => (s[k] === l ? 1 : 0)))); }
  });
  return samples.map((_, i) => cols.map(c => c[i]));
}
function sampleCorrelation(M, n, k) {
  const C = Array.from({ length: n }, () => new Array(n).fill(1));
  const rows = Array.from({ length: n }, (_, i) => M.subarray(i * k, i * k + k));
  for (let i = 0; i < n; i++) for (let l = i + 1; l < n; l++) { const r = E.pearson(rows[i], rows[l]).r; C[i][l] = C[l][i] = r; }
  return C;
}
function sampleDistances(Y, n, k) {
  const out = []; for (let i = 0; i < n; i++) for (let l = i + 1; l < n; l++) { let s = 0; for (let j = 0; j < k; j++) { const d = Y[i * k + j] - Y[l * k + j]; s += d * d; } out.push(Math.sqrt(s)); }
  return out;
}
function variancePartition(p, samples, g, pairs, params) {
  const n = samples.length, cols = p.keep.slice().sort((a, b) => p.featureVar[b] - p.featureVar[a]).slice(0, 400);
  if (pairs) { // balanced two-way decomposition on complete pairs: tissue + individual + residual
    const cl = (DS.covariates.find(c => c.key === DS.contrast.key) || {}).label || 'Group';
    const out = { labels: [cl, DS.pairing.label, 'Residual'], values: [[], [], []] };
    const idx = new Map(); samples.forEach((s, i) => { const k = pairs[i]; if (!idx.has(k)) idx.set(k, [-1, -1]); idx.get(k)[g[i]] = i; });
    for (const j of cols) {
      const prs = []; idx.forEach(([a, b]) => { const x = p.X[a * p.d + j], y = p.X[b * p.d + j]; if (a >= 0 && b >= 0 && x === x && y === y) prs.push([x, y]); });
      const m = prs.length; if (m < 5) continue;
      const all = prs.flat(), gm = E.mean(all); let sst = 0; all.forEach(v => sst += (v - gm) ** 2); if (!(sst > 0)) continue;
      const m0 = E.mean(prs.map(x => x[0])), m1 = E.mean(prs.map(x => x[1]));
      const ssT = m * ((m0 - gm) ** 2 + (m1 - gm) ** 2); let ssI = 0; prs.forEach(([x, y]) => ssI += 2 * ((x + y) / 2 - gm) ** 2);
      const ssR = Math.max(0, sst - ssT - ssI);
      out.values[0].push(ssT / sst); out.values[1].push(ssI / sst); out.values[2].push(ssR / sst);
    }
    return out;
  }
  const keys = (params.vpCovariates || DS.vp_covariates || []).filter(k => DS.covariates.find(c => c.key === k));
  const factors = [{ key: DS.contrast.key, label: DS.covariates.find(c => c.key === DS.contrast.key).label, codes: samples.map(s => s[DS.contrast.key]) }]
    .concat(keys.filter(k => k !== DS.contrast.key).map(k => { const c = DS.covariates.find(x => x.key === k); return { key: k, label: c.label, numeric: c.type === 'numeric', codes: samples.map(s => (c.type === 'numeric' ? +s[k] : s[k])) }; }));
  const r = E.variancePartition(p.X, n, p.d, cols, factors);
  return { labels: r.factors.map(f => f.label).concat(['Residual']), values: r.factors.map(f => f.values).concat([r.residual]) };
}

// ---------- K experiment ----------
function kExperiment(Ks, id) {
  if (!LAST) throw new Error('Run the analysis first.');
  const n = LAST.samples.length, fits = [];
  Ks.forEach(K => { progress(id, `Fitting K = ${K}`); const f = E.gfa(LAST.prep.map(p => ({ M: p.M, k: p.k })), n, K, { seed: 557 }); fits.push({ K, r2: f.r2, Z: f.Z, iterations: f.iterations }); });
  const ref = fits[fits.length - 1];
  return fits.map(f => ({ K: f.K, r2: f.r2, iterations: f.iterations, match: E.matchFactors(f.Z, ref.Z).match.map(x => (x ? Math.abs(x.r) : null)) }));
}

// ---------- stability (refit on subsamples, match factors to the full fit) ----------
function stability(q, id) {
  if (!LAST) throw new Error('Run the analysis first.');
  const { prep, samples, fit } = LAST, n = samples.length, B = q.B || 8, frac = q.frac || 0.8, rand = E.rng(77), K = fit.K;
  const res = Array.from({ length: K }, () => []);
  for (let b = 0; b < B; b++) {
    progress(id, `Refitting on subsample ${b + 1} of ${B}`);
    const keep = E.shuffleIdx(n, rand).slice(0, Math.round(n * frac)).sort((a, c) => a - c);
    const views = prep.map(p => ({ M: E.subsetRows(p.M, n, p.k, keep), k: p.k }));
    const f = E.gfa(views, keep.length, K, { seed: 557 + b });
    const Zfull = keep.map(i => fit.Z[i]);
    const m = E.matchFactors(Zfull, f.Z).match; m.forEach((x, k) => res[k].push(x ? Math.abs(x.r) : 0));
  }
  return { perFactor: res, B, frac };
}

// ---------- enrichment ----------
function enrich(q) {
  if (!LAST || !GS) throw new Error('Gene sets are not loaded yet.');
  const cat = GS.catalogs[q.catalog]; if (!cat) throw new Error('Unknown catalogue');
  const sets = {}; cat.ids.forEach((sid, i) => { sets[sid] = { name: cat.names[i], genes: cat.members[i].map(x => GS.genes[x]) }; });
  const m = q.view, v = DS.views[m];
  // per-gene score from factor weights (model features) or from moderated t (all retained features)
  const score = {};
  if (q.source === 'factor') {
    const W = LAST.fit.W[m]; LAST.prep[m].modelCols.forEach((j, r) => { const s = v.symbols[j]; if (!s) return; const w = W[r][q.factor]; if (!(s in score) || Math.abs(w) > Math.abs(score[s])) score[s] = w; });
  } else {
    const p = LAST.prep[m]; const de = E.featureTests(p.X, LAST.samples.length, p.d, p.keep, LAST.g, null, { minPerGroup: 3, pairs: LAST.pairs }).results;
    de.forEach(r => { const s = v.symbols[r.j]; if (!s || !isFinite(r.t)) return; if (!(s in score) || Math.abs(r.t) > Math.abs(score[s])) score[s] = r.t; });
  }
  const genes = Object.keys(score);
  if (q.method === 'rank') return Object.assign({ method: 'rank' }, E.rankTest(score, sets, q.minSize || 5, q.maxSize || 400));
  const sorted = genes.slice().sort((a, b) => (q.sign === 'pos' ? score[b] - score[a] : q.sign === 'neg' ? score[a] - score[b] : Math.abs(score[b]) - Math.abs(score[a])));
  const sel = sorted.slice(0, q.top || 40);
  return Object.assign({ method: 'ora', selected: sel }, E.ora(sel, genes, sets, q.minSize || 3, q.maxSize || 400));
}

// ---------- network ----------
function network(q) {
  if (!LAST) throw new Error('Run the analysis first.');
  const { prep, fit } = LAST, n = LAST.samples.length, k = q.factor, per = q.per || 8;
  const nodes = [{ id: 'F', type: 'factor', label: `F${k + 1}`, view: -1 }], links = [], chosen = [];
  prep.forEach((p, m) => {
    const W = fit.W[m]; const ord = W.map((w, r) => [r, w[k]]).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1])).slice(0, per);
    ord.forEach(([r, w]) => { const j = p.modelCols[r]; const idn = `${m}:${j}`; nodes.push({ id: idn, view: m, j, r, label: DS.views[m].labels[j], w }); chosen.push({ id: idn, m, r }); if (q.mode === 'factor') links.push({ source: 'F', target: idn, value: w, kind: 'factor' }); });
  });
  if (q.mode === 'corr') {
    nodes.shift();
    const col = (m, r) => { const p = prep[m]; const c = new Float64Array(n); for (let i = 0; i < n; i++) c[i] = p.M[i * p.k + r]; return c; };
    const cols = chosen.map(c => col(c.m, c.r)); const cand = [];
    for (let a = 0; a < chosen.length; a++) for (let b = a + 1; b < chosen.length; b++) {
      if (chosen[a].m === chosen[b].m && !q.within) continue;
      const s = E.spearman(cols[a], cols[b]); if (s.n >= 10 && Math.abs(s.r) >= q.minRho) cand.push({ source: chosen[a].id, target: chosen[b].id, value: s.r, n: s.n, p: s.p, kind: 'corr' });
    }
    cand.sort((x, y) => Math.abs(y.value) - Math.abs(x.value)); links.push(...cand.slice(0, 120));
  }
  return { nodes, links };
}

// ---------- prediction ----------
function predict(q, id) {
  if (!LAST) throw new Error('Run the analysis first.');
  const { prep, samples, g } = LAST, n = samples.length;
  const blocks = prep.map((p, m) => ({ name: DS.views[m].short, M: E.subsetCols(p.X, n, p.d, p.keep), k: p.keep.length }));
  // early integration: concatenate the layers
  const kt = blocks.reduce((a, b) => a + b.k, 0), C = new Float64Array(n * kt); let off = 0;
  blocks.forEach(b => { for (let i = 0; i < n; i++) for (let j = 0; j < b.k; j++) C[i * kt + off + j] = b.M[i * b.k + j]; off += b.k; });
  const all = blocks.concat([{ name: 'All layers', M: C, k: kt }]);
  const y = q.target ? targetLabels(q.target, samples) : { y: g, keep: samples.map(() => true) };
  const rowsKeep = y.keep.map((k, i) => (k ? i : -1)).filter(i => i >= 0);
  const sub = all.map(b => ({ name: b.name, M: E.subsetRows(b.M, n, b.k, rowsKeep), k: b.k }));
  const yy = rowsKeep.map(i => y.y[i]);
  const opt = { folds: q.folds || 5, repeats: q.repeats || 5, topN: q.topN || 20, lambda: q.lambda || 1, seed: 557 };
  progress(id, 'Cross-validating: honest pipeline');
  const honest = E.crossValidate(sub, yy, rowsKeep.length, Object.assign({}, opt));
  progress(id, 'Cross-validating: feature selection before splitting');
  const leaky = E.crossValidate(sub, yy, rowsKeep.length, Object.assign({}, opt, { leaky: true }));
  progress(id, 'Cross-validating: shuffled labels');
  const permH = E.crossValidate(sub, yy, rowsKeep.length, Object.assign({}, opt, { permute: true }));
  const permL = E.crossValidate(sub, yy, rowsKeep.length, Object.assign({}, opt, { permute: true, leaky: true }));
  return { honest, leaky, permH, permL, n: rowsKeep.length, nPos: yy.filter(Boolean).length, target: y.label || DS.contrast.key };
}
function targetLabels(key, samples) {
  const c = DS.covariates.find(x => x.key === key); const vals = samples.map(s => s[key]);
  const lv = [...new Set(vals.filter(v => v != null && v !== '' && v !== 'Normal tissue'))].sort();
  const pos = c && c.positive ? c.positive : ['Mutant', 'Ever', 'Group e', 'Tumor'].find(x => lv.includes(x)) || lv[lv.length - 1];
  return { y: vals.map(v => (v === pos ? 1 : 0)), keep: vals.map(v => v != null && v !== '' && v !== 'Normal tissue'), label: `${c ? c.label : key}: ${pos}` };
}

// ---------- RNA-protein concordance on model input values ----------
function concordance(q) {
  if (!LAST) throw new Error('Run the analysis first.');
  const { prep, samples } = LAST, n = samples.length; const a = q.a, b = q.b;
  const va = DS.views[a], vb = DS.views[b];
  const ia = va.symbols.indexOf(q.gene), ib = vb.symbols.indexOf(q.gene);
  if (ia < 0 || ib < 0) return null;
  const pa = prep[a], pb = prep[b];
  const x = Array.from({ length: n }, (_, i) => pa.X[i * pa.d + ia]), y = Array.from({ length: n }, (_, i) => pb.X[i * pb.d + ib]);
  return { x, y, s: E.spearman(x, y) };
}
