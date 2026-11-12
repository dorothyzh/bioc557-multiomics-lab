/* Multi-omics Lab — UI */
(function () {
  'use strict';
  const C = window.Charts, E = window.Engine, d3 = window.d3;
  const DATASETS = [
    { key: 'honeybee', file: 'data/honeybee_hbb.json', gs: 'data/honeybee_hbb_genesets.json', label: 'Honey bee', sub: 'head · abdomen · gut' },
    { key: 'human', file: 'data/human_luad.json', gs: 'data/human_luad_genesets.json', label: 'Human', sub: 'RNA · protein · phospho' }
  ];
  const VIEW_COLORS = ['#7151A8', '#1B8A6B', '#C2831C', '#B8436F'];
  const GROUP_COLORS = ['#0072B2', '#D55E00'];
  const CAT = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#8c8c8c'];
  const SHAPES = ['circle', 'triangle', 'square', 'diamond'];
  const CHAPTERS = [
    { id: 'study', n: '01', title: 'The study', kicker: 'The study' },
    { id: 'quality', n: '02', title: 'Quality & preparation', kicker: 'Quality & preparation' },
    { id: 'layers', n: '03', title: 'Each layer on its own', kicker: 'Single-layer analysis' },
    { id: 'concordance', n: '04', title: 'Across layers, gene by gene', kicker: 'Concordance', isNew: true },
    { id: 'integration', n: '05', title: 'Integration', kicker: 'Integration' },
    { id: 'network', n: '06', title: 'Network explorer', kicker: 'Network' },
    { id: 'pathways', n: '07', title: 'GO & pathways', kicker: 'Functional interpretation' },
    { id: 'predict', n: '08', title: 'Prediction & leakage', kicker: 'Prediction', isNew: true },
    { id: 'challenge', n: '09', title: 'Break it, then explain it', kicker: 'The challenge' },
    { id: 'export', n: '10', title: 'Export & continue', kicker: 'Your analysis' }
  ];
  const DATA_VERSION = '2026-11-11c';
  const DEFAULTS = { batchAdjust: false, sampleNorm: 'none', maxMissing: 0.3, modelFeatures: 200, scaling: 'unit', center: true, K: 6, impute: 'median', paired: true, challenge: 'none', seed: 557 };

  const S = { ds: null, dsKey: null, gs: null, params: Object.assign({}, DEFAULTS), res: null, chapter: 'study', view: 0, colorBy: null, factor: 0, fx: 0, fy: 1,
    feature: null, figs: new Map(), challengeRes: null, baseline: null, kexp: null, stab: null, pred: null, enrichRes: null, net: null };
  const $ = id => document.getElementById(id);
  const esc = C.esc, fmt = C.fmt, fmtP = C.fmtP;
  const pct = v => (v == null || !isFinite(v) ? '—' : (100 * v).toFixed(1) + '%');

  // ---------- worker ----------
  const worker = new Worker('js/worker.js?v=6'); let reqId = 0; const pending = new Map();
  worker.onmessage = ev => {
    const m = ev.data;
    if (m.type === 'progress') { setStatus('busy', m.msg); return; }
    const p = pending.get(m.id); if (!p) return; pending.delete(m.id);
    if (m.type === 'error') p.reject(new Error(m.message)); else p.resolve(m.result);
  };
  const call = (type, payload) => new Promise((resolve, reject) => { const id = ++reqId; pending.set(id, { resolve, reject }); worker.postMessage(Object.assign({ type, id }, payload)); });
  function setStatus(kind, text) { const el = $('engineStatus'); el.className = 'engine ' + kind; el.querySelector('.txt').textContent = text; }

  // ---------- helpers ----------
  function h(tag, attrs, html) { const el = document.createElement(tag); if (attrs) for (const [k, v] of Object.entries(attrs)) { if (k === 'class') el.className = v; else if (k.startsWith('on')) el.addEventListener(k.slice(2), v); else if (v != null) el.setAttribute(k, v); } if (html != null) el.innerHTML = html; return el; }
  const slug = s => String(s).toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '').slice(0, 40);
  const ds = () => S.ds;
  const contrast = () => S.ds.contrast;
  const groupColor = g => GROUP_COLORS[Math.max(0, contrast().levels.indexOf(g))];
  const groupsSpec = () => contrast().levels.map((l, i) => ({ name: l, label: contrast().labels[i], color: GROUP_COLORS[i], shape: SHAPES[i] }));
  const sampleById = id => S.ds.samples.find(s => s.id === id);
  function covariate(key) { return S.ds.covariates.find(c => c.key === key); }
  function colorSpec(key, ids) { // returns {groups, colorOf(id), shapeOf(id), legend}
    const cov = covariate(key);
    if (!cov || key === contrast().key) return { groups: groupsSpec(), get: s => ({ group: s[contrast().key] }) };
    if (cov.type === 'numeric') {
      const vals = ids.map(id => +sampleById(id)[key]).filter(isFinite); const sc = C.sequential(d3.extent(vals), 'blue');
      return { groups: [], numeric: true, scale: sc, domain: d3.extent(vals), get: s => ({ group: '', color: isFinite(+s[key]) ? sc(+s[key]) : '#ccc' }) };
    }
    const lv = [...new Set(ids.map(id => sampleById(id)[key]))].filter(v => v != null).sort();
    const big = lv.length > CAT.length;
    return { groups: lv.map((l, i) => ({ name: l, label: l, color: big ? d3.interpolateRainbow(i / lv.length) : CAT[i % CAT.length], shape: SHAPES[0] })), many: big, get: s => ({ group: s[key] }) };
  }
  let figCounter = 0;
  function panel(parent, o) {
    const letter = String.fromCharCode(97 + figCounter++);
    const el = h('section', { class: 'panel' + (o.wide ? ' wide' : '') });
    const head = h('div', { class: 'panel-head' }, `<span class="panel-letter">${letter}</span><span class="panel-title">${esc(o.title)}</span>`);
    const tools = h('div', { class: 'panel-tools' });
    const name = () => `${S.ds.id}_fig${chapterMeta().n}${letter}_${slug(o.title)}`;
    tools.append(h('button', { title: 'Download vector figure (SVG, physical size in mm)', onclick: () => { const s = body.querySelector('svg'); if (s) C.exportSVG(s, name()); } }, 'SVG'));
    tools.append(h('button', { title: 'Download 600 dpi PNG', onclick: () => { const s = body.querySelector('svg'); if (s) C.exportPNG(s, name(), 600); } }, 'PNG'));
    head.append(tools); el.append(head);
    const body = h('div', { class: 'panel-body' }); el.append(body);
    const cap = h('div', { class: 'caption' }, o.caption || ''); el.append(cap);
    parent.append(el);
    const draw = (fn) => { try { fn(body); const s = body.querySelector('svg'); if (s) S.figs.set(name(), () => C.serialize(s)); } catch (e) { console.error(e); body.innerHTML = `<div class="empty">Could not draw this figure: ${esc(e.message)}</div>`; } };
    return { el, body, cap, letter, draw, setCaption: html => (cap.innerHTML = html) };
  }
  function pause(parent, q, a) { parent.append(h('div', { class: 'pause' }, `<div class="tag">Pause &amp; explain</div><p>${q}</p><details><summary>Compare your reasoning</summary>${a}</details>`)); }
  function seg(options, value, onChange) {
    const el = h('div', { class: 'seg', role: 'group' });
    options.forEach(o => { const b = h('button', { 'aria-pressed': String(o.value === value), onclick: () => { el.querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', 'false')); b.setAttribute('aria-pressed', 'true'); onChange(o.value); } }, esc(o.label)); el.append(b); });
    return el;
  }
  function select(label, options, value, onChange) {
    const l = h('label', null, esc(label)); const s = h('select', { onchange: () => onChange(s.value) });
    options.forEach(o => { const op = h('option', { value: o.value }, esc(o.label)); if (String(o.value) === String(value)) op.selected = true; s.append(op); });
    l.append(s); return l;
  }
  function viewSeg(onChange) { return seg(S.ds.views.map((v, i) => ({ value: i, label: v.short })), S.view, v => { S.view = +v; onChange(); }); }
  function colorBySelect(onChange) {
    const opts = S.ds.covariates.filter(c => c.type !== 'id' && (S.ds.samples.length < 150 || c.key !== 'site')).map(c => ({ value: c.key, label: c.label }));
    return select('Colour samples by', opts, S.colorBy || contrast().key, v => { S.colorBy = v; onChange(); });
  }
  function stats(parent, items) { const el = h('div', { class: 'stats' }); items.forEach(([v, l]) => el.append(h('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`))); parent.append(el); }
  function sec(parent, title, sub) { parent.append(h('h2', { class: 'sec' }, title)); if (sub) parent.append(h('p', { class: 'sec-sub' }, sub)); }
  function chapterMeta() { return CHAPTERS.find(c => c.id === S.chapter); }
  function orderedSampleIdx(res) { // group, then id
    return res.samples.map((id, i) => i).sort((a, b) => (res.codes[a] - res.codes[b]) || (res.samples[a] < res.samples[b] ? -1 : 1));
  }
  function csv(rows, cols) { const q = v => { const s = v == null ? '' : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }; return [cols.map(c => c[0]).join(',')].concat(rows.map(r => cols.map(c => q(typeof c[1] === 'function' ? c[1](r) : r[c[1]])).join(','))).join('\n'); }
  function download(name, text, type) { const a = h('a', { href: URL.createObjectURL(new Blob([text], { type: type || 'text/csv' })), download: name }); document.body.append(a); a.click(); setTimeout(() => a.remove(), 300); }
  const sigQ = 0.05;

  // ---------- shell ----------
  function buildShell() {
    const sw = $('dsSwitch'); sw.innerHTML = '';
    DATASETS.forEach(d => sw.append(h('button', { role: 'tab', 'aria-selected': String(d.key === S.dsKey), onclick: () => loadDataset(d.key) }, `${d.label}<small>${d.sub}</small>`)));
    const toc = $('toc'); toc.innerHTML = '';
    CHAPTERS.forEach(c => toc.append(h('a', { href: '#' + c.id, class: c.id === S.chapter ? 'active' : '' }, `<span class="n">${c.n}</span><span>${esc(c.title)}${c.isNew ? '<span class="new">new</span>' : ''}</span>`)));
    $('crumb').innerHTML = S.ds ? `<b>${esc(S.ds.species_short)}</b> · ${esc(S.ds.title)}` : '';
  }
  function buildSettings() {
    const el = $('settings'); el.innerHTML = ''; const p = S.params;
    const add = (label, key, opts, conv) => el.append(select(label, opts, p[key], v => { p[key] = conv ? conv(v) : v; scheduleRun(); }));
    add('Sample normalisation', 'sampleNorm', [{ value: 'none', label: 'As supplied' }, { value: 'median', label: 'Median-centre each sample' }]);
    add('Max. missing per feature', 'maxMissing', [0.2, 0.3, 0.5].map(v => ({ value: v, label: pct(v).replace('.0', '') })), Number);
    add('Model features per layer', 'modelFeatures', [100, 200, 300, 500].map(v => ({ value: v, label: v })), Number);
    add('Feature scaling', 'scaling', [{ value: 'unit', label: 'Unit variance' }, { value: 'pareto', label: 'Pareto' }, { value: 'none', label: 'None (centre only)' }]);
    add('Factors (K)', 'K', [3, 4, 5, 6, 7, 8].map(v => ({ value: v, label: v })), Number);
    add('Imputation for PCA', 'impute', [{ value: 'median', label: 'Feature median' }, { value: 'min', label: 'Feature minimum (left-censored)' }]);
    if (S.ds.batch) add('Batch adjustment', 'batchAdjust', [{ value: 'false', label: 'None' }, { value: 'true', label: `Centre features within ${(covariate(S.ds.batch) || {}).label || S.ds.batch}` }], v => v === 'true');
    if (S.ds.pairing) add('Differential tests', 'paired', [{ value: 'true', label: `Paired by ${S.ds.pairing.label.toLowerCase()}` }, { value: 'false', label: 'Unpaired' }], v => v === 'true');
    el.append(h('div', { class: 'actions' }, `<button class="ghost" id="resetSettings">Reset to defaults</button><span>Changes re-run the whole pipeline in your browser. Every chapter updates.</span>`));
    $('resetSettings').onclick = () => { S.params = Object.assign({}, DEFAULTS); buildSettings(); scheduleRun(); };
  }
  let runTimer = null;
  function scheduleRun() { clearTimeout(runTimer); runTimer = setTimeout(runAnalysis, 250); }
  async function runAnalysis() {
    const t0 = performance.now(); setStatus('busy', 'Running the analysis');
    try {
      const res = await call('run', { params: Object.assign({}, S.params, { challenge: 'none' }) });
      S.res = res; S.kexp = null; S.stab = null; S.pred = null; S.enrichRes = null; S.net = null; S.challengeRes = null;
      if (!S.baseline) S.baseline = res;
      S.factor = Math.min(S.factor, res.gfa.K - 1); S.fy = Math.min(S.fy, res.gfa.K - 1);
      setStatus('ready', `Analysis ready · ${(performance.now() - t0) / 1000 < 1 ? '<1' : ((performance.now() - t0) / 1000).toFixed(1)} s`);
      const p = S.params; $('runNote').hidden = false;
      $('runNote').textContent = `${res.n} samples · ${S.ds.views.length} layers · ${p.modelFeatures} model features/layer · ≤${pct(p.maxMissing).replace('.0', '')} missing · ${p.scaling} scaling · K = ${res.gfa.K}${S.ds.pairing ? (p.paired ? ' · paired tests' : ' · unpaired tests') : ''}${S.ds.batch ? (p.batchAdjust ? ' · batch-adjusted' : ' · no batch adjustment') : ''}`;
      render();
    } catch (e) { setStatus('error', 'Analysis failed: ' + e.message); console.error(e); }
  }

  async function loadDataset(key) {
    if (S.dsKey === key && S.ds) return;
    const d = DATASETS.find(x => x.key === key); S.dsKey = key; S.ds = null; S.res = null; S.baseline = null; S.gs = null; S.figs.clear();
    try { localStorage.setItem('mol-ds', key); } catch (e) { /* storage unavailable */ }
    buildShell(); setStatus('busy', `Loading ${d.label.toLowerCase()} data`);
    $('chapter').innerHTML = '<div class="empty">Loading data…</div>';
    const ds = await fetch(d.file + '?v=' + DATA_VERSION).then(r => r.json());
    S.ds = ds; S.view = 0; S.colorBy = ds.contrast.key; S.factor = 0; S.fx = 0; S.fy = 1; S.feature = null;
    S.params = Object.assign({}, DEFAULTS, { K: ds.samples.length < 60 ? 4 : 6 });
    buildShell(); buildSettings(); render();
    await call('load', { dataset: ds });
    fetch(d.gs + '?v=' + DATA_VERSION).then(r => r.json()).then(gs => { S.gs = gs; return call('genesets', { genesets: gs }); }).catch(() => {});
    await runAnalysis();
  }

  // ---------- routing ----------
  function route() {
    const id = (location.hash || '#study').slice(1); S.chapter = CHAPTERS.some(c => c.id === id) ? id : 'study';
    document.querySelectorAll('#toc a').forEach(a => a.classList.toggle('active', a.getAttribute('href') === '#' + S.chapter));
    $('side').classList.remove('open'); render(); window.scrollTo({ top: 0 });
  }
  function render() {
    const root = $('chapter'); root.innerHTML = ''; figCounter = 0;
    const meta = chapterMeta(); const i = CHAPTERS.indexOf(meta);
    if (!S.ds) return;
    const lede = LEDES[meta.id](S.ds);
    root.append(h('div', { class: 'kicker' }, `${meta.n} · ${esc(meta.kicker)}`), h('h1', null, lede[0]), h('p', { class: 'lede' }, lede[1]));
    const tools = h('div', { class: 'chapter-tools' });
    tools.append(h('button', { class: 'ghost', onclick: () => showCode(meta) }, 'R code for this chapter'));
    root.append(tools);
    if (!S.res && meta.id !== 'study') { root.append(h('div', { class: 'empty' }, 'Running the analysis…')); return; }
    try { RENDER[meta.id](root); } catch (e) { console.error(e); root.append(h('div', { class: 'empty' }, 'Something went wrong while drawing this chapter: ' + esc(e.message))); }
    const pager = h('div', { class: 'pager' });
    pager.append(i > 0 ? h('a', { class: 'ghost', href: '#' + CHAPTERS[i - 1].id, style: 'text-decoration:none' }, '← ' + esc(CHAPTERS[i - 1].title)) : h('span'));
    pager.append(i < CHAPTERS.length - 1 ? h('a', { class: 'primary', href: '#' + CHAPTERS[i + 1].id, style: 'text-decoration:none' }, esc(CHAPTERS[i + 1].title) + ' →') : h('span'));
    root.append(pager);
  }
  function showCode(meta) { $('codeTitle').textContent = `R · ${meta.title}`; $('codeBody').textContent = (window.RCODE || {})[meta.id] || ''; $('codeDialog').showModal(); }

  const isHuman = () => S.ds.id === 'human_luad';
  const LEDES = {
    study: d => isHuman() ? ['Three molecular layers, one tumour.', `RNA, proteins and phosphosites from ${d.samples.length / 2} lung adenocarcinomas and the normal lung next to each one.`]
      : ['One pool of bees, three tissues.', `Head, abdomen and midgut proteomes from ${d.samples.length} pools of nurse bees at highbush blueberry apiaries, sampled twice in a season.`],
    quality: () => ['Start with the measurements.', 'Check missingness, abundance scales and sample-level structure before anything is combined.'],
    layers: () => ['What does each layer say on its own?', 'Sample structure, moderated two-group tests and where the variance comes from, one layer at a time.'],
    concordance: () => isHuman() ? ['Does more RNA mean more protein?', 'Follow the same gene across layers and ask which genes are coupled, and which are not.']
      : ['Is a protein high in every tissue at once?', 'Follow the same protein across tissues of the same bees and ask which ones move together.'],
    integration: () => ['Put the layers together.', 'Compare an equally weighted joint PCA with a Bayesian multi-view factor model.'],
    network: () => ['Connect factors and features.', 'Learned factor weights and observed co-variation across layers, as an explorable network.'],
    pathways: () => ['From features to biological questions.', 'Rank-based and over-representation tests on factor weights or differential statistics.'],
    predict: () => ['Can the layers predict the group, honestly?', 'Cross-validated classifiers per layer and combined, and what feature selection outside the folds does to the result.'],
    challenge: () => ['How much does pairing matter?', 'Change one part of the experiment and follow its effect through the same pipeline.'],
    export: () => ['Take the analysis with you.', 'Download figures at journal size, result tables, the exact settings and a companion R walkthrough.']
  };

  // ---------- chapters ----------
  const RENDER = {};

  RENDER.study = root => {
    const d = S.ds;
    const prose = h('div', { class: 'prose' });
    prose.innerHTML = isHuman()
      ? `<p>The Clinical Proteomic Tumor Analysis Consortium (CPTAC) profiled lung adenocarcinomas and normal-adjacent lung from the same patients. We use the harmonized 2023 release, in which every cohort was re-processed with one pipeline. Each patient contributes two samples, so <b>tumour versus normal</b> is a <b>paired</b> comparison.</p><p>Three layers are matched by sample ID: gene-level RNA-seq, TMT mass-spectrometry proteomics, and IMAC-enriched phosphoproteomics. Sex, age, smoking history, stage and driver mutations are available as covariates, which lets you ask what a learned factor actually tracks.</p>`
      : `<p>Colonies were kept at ten highbush blueberry apiaries over two seasons (2020 and 2021). Each apiary had two colony sets (coded <i>e</i> and <i>u</i> in the public files), and each set was sampled at two timepoints (<i>t2</i> and <i>t4</i>). Twelve nurse bees per sample were pooled and dissected into <b>head</b>, <b>abdomen</b> and <b>midgut</b>, and each tissue was measured by data-independent acquisition (DIA) mass spectrometry.</p><p>Here the three views are three <b>tissues</b> from the same bees rather than three molecular types, so integration asks whether a colony-level state shows up across the body. The main contrast is <b>t4 versus t2</b>, paired within each colony set. Season coincides with the mass-spectrometry batch (each year was run and searched separately): a built-in lesson in batch effects that you can switch on and off in the analysis settings.</p>`;
    root.append(prose);
    const cards = h('div', { class: 'cards' });
    DATASETS.forEach(x => {
      const active = x.key === S.dsKey;
      const c = h('div', { class: 'ds-card' + (active ? ' active' : ''), onclick: () => loadDataset(x.key) });
      c.innerHTML = x.key === 'human'
        ? `<div class="sp">Homo sapiens · CPTAC 2023</div><h3>Lung adenocarcinoma and matched normal lung</h3><div class="layer-row"><span class="chip"><i style="background:${VIEW_COLORS[0]}"></i>RNA-seq</span><span class="chip"><i style="background:${VIEW_COLORS[1]}"></i>TMT proteome</span><span class="chip"><i style="background:${VIEW_COLORS[2]}"></i>Phosphoproteome</span></div><p>101 patients · 202 samples · paired tumour/normal design · sex, age, smoking, EGFR/KRAS/TP53/STK11 status.</p><div class="cite">Li et al., Cancer Cell 2023; Gillette et al., Cell 2020</div>`
        : `<div class="sp">Apis mellifera · 2025</div><h3>Nurse bees at highbush blueberry apiaries</h3><div class="layer-row"><span class="chip"><i style="background:${VIEW_COLORS[0]}"></i>Head proteome</span><span class="chip"><i style="background:${VIEW_COLORS[1]}"></i>Abdomen proteome</span><span class="chip"><i style="background:${VIEW_COLORS[2]}"></i>Gut proteome</span></div><p>40 pools of 12 bees · 10 apiaries × 2 colony sets × 2 timepoints · two seasons run as separate MS batches · Nosema and virus protein signal as covariates.</p><div class="cite">Zhong et al., Proteomics 2025 · MassIVE MSV000097590</div>`;
      cards.append(c);
    });
    root.append(cards);
    if (!S.res) return;
    stats(root, [[d.samples.length, 'matched samples'], [d.views.length, 'layers'], ...d.views.map((v, i) => [S.res.views[i].qc.dInput.toLocaleString(), `${esc(v.short)} features shipped`])]);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    // a: coverage heatmap (observed fraction per sample and layer)
    const ord = orderedSampleIdx(S.res);
    const a = panel(grid, { title: 'Measurement coverage per sample and layer', wide: true, caption: `Fraction of shipped features observed in each sample (columns, ordered by ${esc(contrast().key)}). Coloured strip: ${groupsSpec().map(g => `<b style="color:${g.color}">${esc(g.label)}</b>`).join(' / ')}. Blocks of lower coverage in mass-spectrometry layers usually trace back to multiplexed batches or runs.` });
    a.draw(body => C.heatmap(body, { width: 'double', rows: S.res.views.map(v => v.short), cols: ord.map(i => S.res.samples[i]), values: S.res.views.map(v => ord.map(i => 1 - v.qc.sampleMissing[i])),
      scale: 'sequential', domain: [0.5, 1], hue: 'teal', cellH: 14, ml: 52, mb: 10, showColLabels: false, colAnnot: ord.map(i => ({ color: GROUP_COLORS[S.res.codes[i]] })), legendTitle: 'Observed',
      tooltip: c => `<b>${esc(S.res.samples[ord[c.j]])}</b><br>${esc(S.res.views[c.i].name)}: ${pct(c.v)} observed` }));
    // b: design table
    const second = isHuman() ? 'sex' : 'group', third = isHuman() ? 'smoking' : 'year';
    const lv2 = [...new Set(d.samples.map(s => s[second]))].filter(Boolean).sort(), lv3 = [...new Set(d.samples.map(s => s[third]))].filter(Boolean).sort();
    const cols = []; lv2.forEach(a2 => lv3.forEach(a3 => cols.push([a2, a3])));
    const b = panel(grid, { title: 'Study design', caption: `Samples per ${esc(contrast().key)} level and combination of ${esc(covariate(second).label.toLowerCase())} and ${esc(covariate(third).label.toLowerCase())}.` });
    b.draw(body => C.heatmap(body, { rows: contrast().labels, cols: cols.map(c => c.join(' · ')), values: contrast().levels.map(l => cols.map(([a2, a3]) => d.samples.filter(s => s[contrast().key] === l && s[second] === a2 && s[third] === a3).length)),
      scale: 'sequential', domain: [0, d3.max(cols.map(([a2, a3]) => d.samples.filter(s => s[second] === a2 && s[third] === a3).length)) || 1], cellLabels: v => v, cell: 30, cellH: 20, ml: 62, mb: 46, legendTitle: 'Samples' }));
    const c = panel(grid, { title: 'What each layer measures', caption: 'Assay, unit and filtering of the shipped teaching matrices. Full provenance is on the Sources page.' });
    c.body.innerHTML = `<table class="data"><thead><tr><th>Layer</th><th>Assay</th><th>Unit</th></tr></thead><tbody>${d.views.map((v, i) => `<tr><td><span class="chip"><i style="background:${VIEW_COLORS[i]}"></i>${esc(v.name)}</span></td><td>${esc(v.assay)}</td><td>${esc(v.unit)}</td></tr>`).join('')}</tbody></table><p class="note">${d.views.map(v => esc(v.note)).join(' ')}</p>`;
    sec(root, 'Read before you analyse', 'Every dataset has limits. These shape what you are allowed to conclude.');
    root.append(h('ul', { class: 'prose' }, d.source.limitations.map(l => `<li>${esc(l)}</li>`).join('')));
  };

  RENDER.quality = root => {
    const v = S.res.views[S.view], dv = S.ds.views[S.view];
    const ctr = h('div', { class: 'controls' }); ctr.append(viewSeg(() => render()), colorBySelect(() => render())); root.append(ctr);
    stats(root, [[v.qc.dInput.toLocaleString(), 'features shipped'], [v.qc.dPass.toLocaleString(), `pass ≤ ${pct(S.params.maxMissing).replace('.0', '')} missing`], [pct(v.qc.missingTotal), 'of measurements missing'], [v.qc.nImputedForPCA.toLocaleString(), 'model entries imputed (PCA only)']]);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const ord = orderedSampleIdx(S.res);
    const legendItems = groupsSpec().map(g => ({ label: g.label, color: g.color }));
    const a = panel(grid, { title: 'Missing values per sample', caption: `Before feature filtering. Samples ordered by group. A sample far above the rest deserves a look at its run before it is removed.` });
    a.draw(body => C.sampleBars(body, { items: ord.map(i => ({ id: S.res.samples[i], value: v.qc.sampleMissing[i], color: GROUP_COLORS[S.res.codes[i]] })), yLabel: 'Fraction missing', xLabel: 'Sample', yFormat: '.0%', legendItems, minMax: 0.05,
      tooltip: d => `<b>${esc(d.id)}</b><br>${pct(d.value)} missing` }));
    const b = panel(grid, { title: 'Is missingness related to abundance?', caption: `Each point is a feature: its mean observed log₂ abundance against the fraction of samples where it is missing. A falling trend means low-abundance features drop out: values are <b>missing not at random</b> (left-censored), typical of mass spectrometry.` });
    b.draw(body => {
      const pts = v.qc.featureMean.map((m, j) => ({ x: m, y: v.qc.featureMissing[j], label: dv.labels[j], group: 'f' })).filter(p => isFinite(p.x));
      const bins = d3.bin().value(p => p.x).thresholds(12)(pts).filter(bb => bb.length >= 5);
      const f = C.scatter(body, { points: pts, groups: [{ name: 'f', color: VIEW_COLORS[S.view] }], xLabel: `Mean ${dv.unit}`, yLabel: 'Fraction missing', r: 1.1, opacity: 0.35, legend: false,
        tooltip: p => `<b>${esc(p.label)}</b><br>mean ${fmt(p.x)} · ${pct(p.y)} missing` });
      const x = d3.scaleLinear().domain(f.svg.select('g').node() ? d3.extent(pts, p => p.x) : [0, 1]);
      // binned trend line drawn on top in plot coordinates
      const xs = d3.scaleLinear().domain(pad(d3.extent(pts, p => p.x))).range([0, f.iw]).nice(), ys = d3.scaleLinear().domain(pad(d3.extent(pts, p => p.y))).range([f.ih, 0]).nice();
      f.g.append('path').attr('d', d3.line()(bins.map(bb => [xs(d3.mean(bb, p => p.x)), ys(d3.mean(bb, p => p.y))]))).attr('fill', 'none').attr('stroke', '#16191d').attr('stroke-width', 1.2);
      return x;
    });
    const c = panel(grid, { title: 'Abundance distribution per sample', caption: `Box: interquartile range; whiskers: 5th–95th percentiles; white tick: median. ${S.params.sampleNorm === 'median' ? 'Samples are median-centred.' : 'Values as supplied by the study.'} Unit: ${esc(dv.unit)}.` });
    c.draw(body => C.sampleBoxes(body, { items: ord.map(i => ({ id: S.res.samples[i], q: v.qc.sampleBox[i], color: GROUP_COLORS[S.res.codes[i]] })), yLabel: dv.unit.length > 26 ? 'log₂ abundance' : dv.unit, xLabel: 'Sample', legendItems,
      tooltip: d => `<b>${esc(d.id)}</b><br>median ${fmt(d.q[2])}<br>IQR ${fmt(d.q[1])} to ${fmt(d.q[3])}` }));
    const dd = panel(grid, { title: 'Sample–sample correlation', caption: `Pearson correlation between samples on the ${S.params.modelFeatures} model features (centred${S.params.scaling === 'unit' ? ', unit variance' : ''}). Ordered by group; the coloured strip marks it. Look for blocks that do not follow the biology.` });
    dd.draw(body => { const m = Math.max(0.2, d3.max(ord.flatMap(i => ord.filter(j => j !== i).map(j => Math.abs(v.corr[i][j] || 0)))) || 1);
      C.heatmap(body, { rows: ord.map(i => S.res.samples[i]), cols: ord.map(i => S.res.samples[i]), values: ord.map(i => ord.map(j => (i === j ? null : v.corr[i][j]))), scale: 'diverging', max: +m.toFixed(2), showRowLabels: false, showColLabels: false, ml: 12, mb: 12, mt: 4,
        colAnnot: ord.map(i => ({ color: GROUP_COLORS[S.res.codes[i]] })), rowAnnot: ord.map(i => ({ color: GROUP_COLORS[S.res.codes[i]] })), legendTitle: 'r', tooltip: cc => `${esc(S.res.samples[ord[cc.i]])} × ${esc(S.res.samples[ord[cc.j]])}<br>r = ${fmt(cc.v)}` }); });
    const e = panel(grid, { title: 'Sample structure in this layer', wide: false, caption: `PCA on the model features; missing entries filled by the feature ${S.params.impute === 'min' ? 'minimum' : 'median'} for PCA only. Colour shows ${esc((covariate(S.colorBy) || {}).label || 'group')}, which does not train the PCA.` });
    e.draw(body => pcaPlot(body, v.pca, 0, 1));
    pause(root, 'Does a sample with many missing values also look unusual in the PCA? What would you check before removing it?',
      '<p>Missingness and PCA position are only loosely linked once missing entries are imputed with a central value: an imputed sample is pulled toward the centre, not outward. Before removing a sample, check its run (identifications, total intensity, chromatography), whether its batch-mates look similar, and whether it is an extreme of a real biological gradient. Removing samples because they "look odd" without a technical reason biases the study.</p>');
    if (S.ds.views.some(x => /proteom|phospho/i.test(x.name))) pause(root, 'Switch the PCA imputation to "feature minimum". Which samples move, and why?', '<p>Filling with the minimum treats every missing value as "below detection". Samples with many missing values are pushed toward the low end of many features at once, so they can separate along a PC that tracks <i>missingness</i> rather than biology. Median filling does the opposite and hides real low-abundance states. Neither is neutral; report which you used.</p>');
  };
  function pad(ext) { const s = (ext[1] - ext[0]) || 1; return [ext[0] - s * 0.06, ext[1] + s * 0.06]; }
  function pcaPlot(body, pc, a, b, extra) {
    const key = S.colorBy || contrast().key, spec = colorSpec(key, S.res.samples);
    const pts = S.res.samples.map((id, i) => { const s = sampleById(id); const g = spec.get(s); return { id, x: pc.scores[a][i], y: pc.scores[b][i], group: g.group, color: g.color, shape: key === contrast().key ? SHAPES[S.res.codes[i]] : undefined, label: id }; });
    return C.scatter(body, Object.assign({ points: pts, groups: spec.groups, xLabel: `PC${a + 1} (${pct(pc.varExp[a])})`, yLabel: `PC${b + 1} (${pct(pc.varExp[b])})`, ellipses: key === contrast().key, zeroLines: true, legend: !spec.many,
      tooltip: p => { const s = sampleById(p.id); return `<b>${esc(p.id)}</b><br>${S.ds.covariates.filter(c => c.type !== 'id').slice(0, 5).map(c => `${esc(c.label)}: ${esc(s[c.key])}`).join('<br>')}`; } }, extra || {}));
  }

  RENDER.layers = root => {
    const v = S.res.views[S.view], dv = S.ds.views[S.view], paired = S.res.paired;
    const ctr = h('div', { class: 'controls' }); ctr.append(viewSeg(() => { S.feature = null; render(); }), colorBySelect(() => render())); root.append(ctr);
    const tested = v.de.filter(r => isFinite(r.p)), nsig = tested.filter(r => r.q < sigQ).length, nsigW = tested.filter(r => r.welchQ < sigQ).length;
    stats(root, [[v.qc.dPass.toLocaleString(), 'features tested'], [nsig.toLocaleString(), `BH q < 0.05 · moderated${paired ? ', paired' : ''}`], [nsigW.toLocaleString(), 'BH q < 0.05 · Welch, unpaired'], [isFinite(v.d0) ? fmt(v.d0, 1) : '∞', 'prior degrees of freedom (eBayes)']]);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const a = panel(grid, { title: `${dv.name}: samples in the first two components`, caption: `PCA on ${S.params.modelFeatures} label-free, high-variance features. Ellipses: 95% normal-theory regions per group.` });
    a.draw(body => pcaPlot(body, v.pca, 0, 1));
    const b = panel(grid, { title: 'Variance captured by each component', caption: 'Share of total variance of the model matrix per principal component. A dominant PC1 means one axis of variation organises this layer.' });
    b.draw(body => C.barsH(body, { items: v.pca.varExp.map((x, i) => ({ label: `PC${i + 1}`, value: x, color: VIEW_COLORS[S.view] })), xLabel: 'Variance explained', xFormat: '.0%', valueFormat: x => pct(x), labelWidth: 30, rowH: 13 }));
    const lv = contrast().labels;
    const sorted = tested.slice().sort((x, y) => x.p - y.p);
    if (!S.feature || S.feature.view !== S.view) S.feature = { view: S.view, j: sorted.length ? sorted[0].j : 0 };
    const c = panel(grid, { title: `Differences between ${lv[1].toLowerCase()} and ${lv[0].toLowerCase()}`, caption: `Moderated t-test (limma-style empirical Bayes${paired ? `, paired by ${esc(S.ds.pairing.label.toLowerCase())}` : ''}) on observed values only; no imputed value enters a test. Coloured: BH q &lt; 0.05. Click a point to inspect that feature.` });
    c.draw(body => {
      const pts = tested.map(r => ({ x: r.diff, y: -Math.log10(Math.max(r.p, 1e-300)), sig: r.q < sigQ, label: r.label, r }));
      const lab = pts.filter(p => p.sig).sort((x, y) => (y.y + 4 * Math.abs(y.x)) - (x.y + 4 * Math.abs(x.x))).slice(0, 8);
      C.volcano(body, { points: pts, labels: lab, xLabel: `Difference (${lv[1]} − ${lv[0]}), ${dv.unit.replace(/ \(.*\)/, '')}`, upColor: GROUP_COLORS[1], downColor: GROUP_COLORS[0], upLabel: lv[1].toLowerCase(), downLabel: lv[0].toLowerCase(), italic: isHuman() && S.view < 2,
        pLine: (() => { const s = tested.filter(r => r.q < sigQ); return s.length ? -Math.log10(d3.max(s, r => r.p)) : null; })(),
        tooltip: p => `<b>${esc(p.label)}</b><br>Δ = ${fmt(p.r.diff)}<br>p = ${fmtP(p.r.p)} · q = ${fmtP(p.r.q)}<br>${paired ? `${p.r.nPairs} complete pairs` : `n = ${p.r.n0} / ${p.r.n1}`}`,
        onClick: p => { S.feature = { view: S.view, j: p.r.j }; drawFeature(); } });
    });
    const dpanel = panel(grid, { title: 'One feature, every sample', caption: '' });
    function drawFeature() {
      const j = S.feature.j, r = v.de.find(x => x.j === j), lab = dv.labels[j];
      dpanel.el.querySelector('.panel-title').textContent = `${lab}: values by group`;
      const ids = S.res.samples, vals = ids.map(id => { const row = dv.sample_ids.indexOf(id); return dv.values[row][j]; });
      const pts = ids.map((id, i) => ({ cat: S.res.groups[i], y: vals[i] == null ? NaN : vals[i], label: id, pair: S.ds.pairing ? sampleById(id)[S.ds.pairing.key] : null })).filter(p => isFinite(p.y));
      dpanel.draw(body => {
        const f = C.strip(body, { categories: contrast().levels.map((l, i) => ({ name: l, label: contrast().labels[i], color: GROUP_COLORS[i] })), points: pts, yLabel: dv.unit.length > 26 ? 'log₂ abundance' : dv.unit,
          stat: r && isFinite(r.p) ? `Δ = ${fmt(r.diff)} · q = ${fmtP(r.q)}` : '' });
        if (paired) { // connect pairs
          const xb = f.g.selectAll('circle').nodes(); const byPair = new Map();
          xb.forEach(nd => { const p = d3.select(nd).datum(); if (!byPair.has(p.pair)) byPair.set(p.pair, []); byPair.get(p.pair).push([+nd.getAttribute('cx'), +nd.getAttribute('cy')]); });
          const lg = f.g.insert('g', 'circle'); byPair.forEach(seg2 => { if (seg2.length === 2) lg.append('line').attr('x1', seg2[0][0]).attr('y1', seg2[0][1]).attr('x2', seg2[1][0]).attr('y2', seg2[1][1]).attr('stroke', '#9aa3ab').attr('stroke-width', 0.3).attr('stroke-opacity', 0.6); });
        }
      });
      dpanel.setCaption(`Observed values of <b>${esc(lab)}</b>${paired ? '; grey lines join the two samples of each pair' : ''}. Box: median and interquartile range. Click any point in the volcano plot or a table row to change the feature.`);
      document.querySelectorAll('#deTable tr').forEach(tr => tr.classList.toggle('sel', +tr.dataset.j === j));
    }
    drawFeature();
    const vp = S.res.vp[S.view];
    const e = panel(grid, { title: 'Where does the variance come from?', wide: false, caption: paired ? `For each of the 400 most variable features, the share of variance between ${esc(contrast().key)} levels, between ${esc(S.ds.pairing.label.toLowerCase())}s (individual differences), and the remainder. Computed on complete pairs.` : `Type-II share of variance for each covariate in a linear model per feature (400 most variable features); the remainder is residual.` });
    const VPC = ['#009E73', '#CC79A7', '#E69F00', '#56B4E9'];
    e.draw(body => C.strip(body, { categories: vp.labels.map((l, i) => ({ name: l, label: l, color: i === vp.labels.length - 1 ? '#8c8c8c' : VPC[i % VPC.length] })), points: vp.labels.flatMap((l, i) => vp.values[i].map(x => ({ cat: l, y: x, label: l }))), yLabel: 'Share of variance', r: 1.1,
      tooltip: p => `${esc(p.cat)}: ${pct(p.y)}` }));
    const t = panel(grid, { title: 'Top features', caption: `Sorted by p. Δ = mean difference (${esc(lv[1])} − ${esc(lv[0])}${paired ? ', within pairs' : ''}). Welch p ignores pairing and does not borrow information across features.` });
    t.body.innerHTML = `<div class="table-wrap"><table class="data" id="deTable"><thead><tr><th>Feature</th><th class="num">Δ</th><th class="num">p</th><th class="num">q</th><th class="num">Welch p</th></tr></thead><tbody>${sorted.slice(0, 60).map(r => `<tr data-j="${r.j}" style="cursor:pointer"><td>${esc(r.label)}</td><td class="num">${fmt(r.diff)}</td><td class="num">${fmtP(r.p)}</td><td class="num">${fmtP(r.q)}</td><td class="num">${fmtP(r.welchP)}</td></tr>`).join('')}</tbody></table></div>`;
    t.body.querySelectorAll('tr[data-j]').forEach(tr => tr.onclick = () => { S.feature = { view: S.view, j: +tr.dataset.j }; drawFeature(); });
    t.el.querySelector('.panel-tools').prepend(h('button', { onclick: () => download(`${S.ds.id}_${dv.key}_differential.csv`, csv(sorted, [['feature', 'f'], ['label', 'label'], ['difference', 'diff'], ['moderated_t', 't'], ['p', 'p'], ['q_BH', 'q'], ['welch_p', 'welchP'], ['welch_q', 'welchQ'], ['n_ref', 'n0'], ['n_test', 'n1'], ['n_pairs', 'nPairs']])) }, 'CSV'));
    drawFeature();
    if (S.ds.pairing) pause(root, `Switch “Differential tests” in the analysis settings between paired and unpaired. Why does the number of discoveries change so much?`, `<p>In a paired design each ${esc(S.ds.pairing.label.toLowerCase())} acts as its own control. The pairing removes between-${esc(S.ds.pairing.label.toLowerCase())} variation (the "Individual" share in the variance plot) from the error term, so the same mean difference is measured against a much smaller noise. Ignoring the pairing does not make a test conservative in a useful way: it answers a different question with less power.</p>`);
    pause(root, 'Do the layers tell the same story? Switch layers and pick one feature you would follow up.', '<p>Expect broad agreement on the strongest contrast but different details: some changes are visible only in one layer (e.g. a phosphosite that moves while its protein does not, or a tissue-specific response). A good follow-up candidate is large, consistent within pairs, measured without many missing values, and plausible in more than one layer.</p>');
  };

  RENDER.concordance = root => {
    const views = S.ds.views; let A = 0, B = 1;
    if (S.concPair) [A, B] = S.concPair;
    const ctr = h('div', { class: 'controls' });
    const pairs = []; for (let i = 0; i < views.length; i++) for (let j = i + 1; j < views.length; j++) pairs.push([i, j]);
    ctr.append(select('Compare', pairs.map(([i, j]) => ({ value: `${i}-${j}`, label: `${views[i].short} ↔ ${views[j].short}` })), `${A}-${B}`, v => { S.concPair = v.split('-').map(Number); render(); }));
    root.append(ctr);
    // per-feature correlations on the shipped pools (browser)
    const va = views[A], vb = views[B], ids = S.res.samples;
    const rowA = new Map(va.sample_ids.map((id, i) => [id, i])), rowB = new Map(vb.sample_ids.map((id, i) => [id, i]));
    const symB = new Map(); vb.symbols.forEach((s, j) => { if (s && !symB.has(s)) symB.set(s, j); });
    const shared = []; va.symbols.forEach((s, j) => { if (s && symB.has(s)) shared.push({ sym: s, ja: j, jb: symB.get(s), label: va.labels[j] }); });
    const vec = (v, rows, j) => ids.map(id => { const x = v.values[rows.get(id)][j]; return x == null ? NaN : x; });
    shared.forEach(g => { const s = E.spearman(vec(va, rowA, g.ja), vec(vb, rowB, g.jb)); g.rho = s.r; g.n = s.n; g.p = s.p; });
    const ok = shared.filter(g => isFinite(g.rho) && g.n >= 10);
    const conc = S.ds.concordance && A === 0 && B === 1 ? S.ds.concordance : null;
    const medShip = d3.median(ok, g => g.rho);
    stats(root, [[ok.length.toLocaleString(), 'features in both teaching pools'], [fmt(medShip), 'median Spearman ρ (pools)'], ...(conc ? [[conc.gene.length.toLocaleString(), 'genes in the full matrices'], [fmt(d3.median(conc.rho_tumor)), 'median ρ, tumours only (full)']] : [])]);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const a = panel(grid, { title: conc ? 'RNA–protein correlation across all measured genes' : `${va.short}–${vb.short} correlation per feature`, caption: conc ? `Per-gene Spearman ρ across the 101 tumours, computed offline on all ${conc.gene.length.toLocaleString()} genes quantified in both layers (≥ 20 pairs). Dashed: median. Within tissue type, so the tumour–normal difference does not inflate ρ.` : `Spearman ρ across all ${ids.length} samples for every feature present in both teaching pools. Dashed: median.` });
    a.draw(body => { const vals = conc ? conc.rho_tumor : ok.map(g => g.rho); const md = d3.median(vals); C.histogram(body, { values: vals, domain: [-1, 1], bins: 40, xLabel: 'Spearman ρ', yLabel: conc ? 'Genes' : 'Features', unit: conc ? 'genes' : 'features', color: VIEW_COLORS[B], markers: [{ value: md, label: `median ${fmt(md)}`, dash: '3,2' }] }); });
    if (conc) {
      const b = panel(grid, { title: 'Abundant proteins track their RNA more closely', caption: 'Genes binned into deciles of mean protein abundance; line: median ρ per bin, band: interquartile range. Low-abundance proteins are measured with more noise and more missing values, which lowers their observable correlation.' });
      b.draw(body => {
        const idx = conc.gene.map((g, i) => i).filter(i => isFinite(conc.protein_mean[i]) && isFinite(conc.rho_tumor[i])).sort((x, y) => conc.protein_mean[x] - conc.protein_mean[y]);
        const k = 10, pts = [], band = []; for (let q = 0; q < k; q++) { const sl = idx.slice(Math.floor(q * idx.length / k), Math.floor((q + 1) * idx.length / k)); const r = sl.map(i => conc.rho_tumor[i]).sort((x, y) => x - y); const xm = d3.median(sl, i => conc.protein_mean[i]); pts.push([xm, d3.quantile(r, 0.5)]); band.push([xm, d3.quantile(r, 0.25), d3.quantile(r, 0.75)]); }
        C.lines(body, { series: [{ label: 'Median ρ', color: VIEW_COLORS[1], points: pts, band, dots: true }], xLabel: 'Mean protein abundance (log₂, decile median)', yLabel: 'Spearman ρ', legend: false, yInclude: [0].concat(band.flatMap(b2 => [b2[1], b2[2]])) });
      });
    }
    // gene picker + scatter
    const strong = ok.slice().sort((x, y) => y.rho - x.rho);
    if (!S.concGene || !ok.find(g => g.sym === S.concGene)) S.concGene = (strong[Math.floor(strong.length * 0.1)] || strong[0] || {}).sym;
    const c = panel(grid, { title: 'Follow one feature across layers', caption: '' });
    const pick = h('div', { class: 'controls' });
    const dl = h('datalist', { id: 'concList' }); ok.slice().sort((x, y) => (x.label < y.label ? -1 : 1)).forEach(g => dl.append(h('option', { value: g.label })));
    const inp = h('input', { type: 'text', list: 'concList', value: (ok.find(g => g.sym === S.concGene) || {}).label || '', placeholder: 'Type a gene or protein', onchange: () => { const g = ok.find(x => x.label === inp.value); if (g) { S.concGene = g.sym; drawGene(); } } });
    const lab = h('label', null, 'Feature'); lab.append(inp); pick.append(lab, dl);
    const quick = h('div', { class: 'seg' }); [['Most coupled', strong[0]], ['Typical', strong[Math.floor(strong.length / 2)]], ['Least coupled', strong[strong.length - 1]]].forEach(([l, g]) => { if (g) quick.append(h('button', { onclick: () => { S.concGene = g.sym; inp.value = g.label; drawGene(); } }, l)); });
    pick.append(quick); c.el.insertBefore(pick, c.body);
    function drawGene() {
      const g = ok.find(x => x.sym === S.concGene); if (!g) return;
      const xa = vec(va, rowA, g.ja), xb = vec(vb, rowB, g.jb);
      const pts = ids.map((id, i) => ({ id, x: xa[i], y: xb[i], group: S.res.groups[i], shape: SHAPES[S.res.codes[i]], label: id }));
      const per = contrast().levels.map(l => { const sel = pts.filter(p => p.group === l); return `${l}: ρ = ${fmt(E.spearman(sel.map(p => p.x), sel.map(p => p.y)).r)}`; }).join(' · ');
      c.draw(body => C.scatter(body, { points: pts, groups: groupsSpec(), xLabel: `${va.short} (${va.unit.length > 22 ? 'log₂' : va.unit})`, yLabel: `${vb.short} (${vb.unit.length > 22 ? 'log₂' : vb.unit})`, fit: true, stat: `all: ρ = ${fmt(g.rho)} · ${per}` }));
      c.setCaption(`<b>${esc(g.label)}</b> in ${esc(va.short)} and ${esc(vb.short)}, one point per sample. Pooled correlations can be driven by the group difference alone; compare with the within-group values above the plot.`);
    }
    drawGene();
    // gene-set test on rho
    if (S.gs) {
      const sc = {}; if (conc) conc.gene.forEach((gname, i) => { if (isFinite(conc.rho_tumor[i])) sc[gname] = conc.rho_tumor[i]; }); else ok.forEach(g => { sc[g.sym] = g.rho; });
      const catKey = S.gs.catalogs.hallmark ? 'hallmark' : 'gobp', cat = S.gs.catalogs[catKey];
      const sets = {}; cat.ids.forEach((id, i) => { sets[id] = { name: cat.names[i], genes: cat.members[i].map(x => S.gs.genes[x]) }; });
      const rt = E.rankTest(sc, sets, 10, 400);
      const top = rt.results.filter(x => x.z > 0).slice(0, 6).concat(rt.results.filter(x => x.z < 0).slice(0, 6));
      const e = panel(grid, { title: 'Which processes are coupled, and which are decoupled?', caption: `Rank test of ρ within ${esc(S.gs.labels[catKey])} sets (Mann–Whitney, BH-adjusted). Effect = AUC − 0.5: positive sets have higher ρ than other genes. Complexes assembled from many subunits (ribosome, proteasome, respiratory chain) are classic examples of post-transcriptional control.` });
      e.draw(body => C.barsH(body, { items: top.sort((x, y) => y.auc - x.auc).map(x => ({ label: x.name.length > 34 ? x.name.slice(0, 33) + '…' : x.name, value: x.auc - 0.5, color: x.z > 0 ? VIEW_COLORS[1] : '#8c8c8c', q: x.q, n: x.size })), domain: [-0.5, 0.5], xLabel: 'Effect on ρ (AUC − 0.5)', labelWidth: 120, rowH: 11, valueFormat: v => fmt(v), tooltip: d => `<b>${esc(d.label)}</b><br>${d.n} genes · q = ${fmtP(d.q)}` }));
    }
    pause(root, isHuman() ? 'The median RNA–protein correlation is around 0.5. Is the missing half "noise"?' : 'Why might a protein correlate between head and abdomen of the same pool?',
      isHuman() ? '<p>Part of it is measurement noise, especially for low-abundance proteins. But a large part is biology: translation rates, protein half-lives, complex stoichiometry (excess subunits are degraded) and secretion all decouple protein from mRNA. That is why proteins carry information that RNA alone does not, and why multi-omics integration is worth doing.</p>'
        : '<p>A shared colony-level state (diet, exposure, infection, age structure of the pool) can raise or lower the same protein throughout the body; secreted or haemolymph proteins such as vitellogenin also move between tissues. Low correlation is expected for tissue-specific enzymes. Watch for batch: year is confounded with the mass-spectrometry run, so a shared "year" shift alone can create correlation across tissues.</p>');
  };

  RENDER.integration = root => {
    const r = S.res, G = r.gfa, K = G.K;
    root.append(h('div', { class: 'prose' }, `<p>The factor model is a Bayesian group factor analysis with automatic relevance determination per layer: the same Gaussian model as MOFA, without the spike-and-slab sparsity prior. It runs in your browser; on the same inputs its factors match MOFA2 (R) with |r| ≥ 0.999. Each layer keeps its own features and noise level. Factors are shared axes of sample variation, and a factor can be active in one layer, several, or all.</p>`));
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const ctr = h('div', { class: 'controls' }); ctr.append(colorBySelect(() => render())); root.insertBefore(ctr, grid);
    const a = panel(grid, { title: 'Joint PCA: the simple baseline', caption: 'All layers concatenated after scaling, each divided by √(number of features) so that no layer dominates by size. Missing values filled with feature medians.' });
    a.draw(body => pcaPlot(body, r.joint, 0, 1));
    const b = panel(grid, { title: 'The factor model: a shared coordinate system', caption: 'Factor scores per sample. Colour and shape show covariates that were not used to fit the model. Factor order and sign are arbitrary.' });
    const fsel = h('div', { class: 'controls' });
    fsel.append(select('x axis', d3.range(K).map(k => ({ value: k, label: `Factor ${k + 1}` })), S.fx, v => { S.fx = +v; drawB(); }), select('y axis', d3.range(K).map(k => ({ value: k, label: `Factor ${k + 1}` })), S.fy, v => { S.fy = +v; drawB(); }));
    b.el.insertBefore(fsel, b.body);
    function drawB() { b.draw(body => pcaPlot(body, { scores: d3.range(K).map(k => G.Z.map(z => z[k])), varExp: d3.range(K).map(() => NaN) }, S.fx, S.fy, { xLabel: `Factor ${S.fx + 1}`, yLabel: `Factor ${S.fy + 1}` })); }
    drawB();
    const c = panel(grid, { title: 'Variance explained per layer and factor', caption: 'R² of each factor alone, computed on observed entries of the centred model matrices. Contributions of different factors can overlap, so they need not add up to the total.' });
    c.draw(body => C.heatmap(body, { rows: r.views.map(v => v.short), cols: d3.range(K).map(k => `F${k + 1}`), values: G.r2.map(x => x.perFactor.map(v => 100 * v)), scale: 'sequential', domain: [0, Math.max(10, d3.max(G.r2.flatMap(x => x.perFactor)) * 100)], cellLabels: v => (v >= 0.95 ? v.toFixed(0) : v.toFixed(1)), cell: 26, cellH: 20, ml: 56, mb: 26, colAngle: 0, legendTitle: 'R² (%)',
      tooltip: cc => `${esc(r.views[cc.i].name)} · Factor ${cc.j + 1}<br>R² = ${fmt(cc.v, 1)}%` }));
    const dd = panel(grid, { title: 'Total variance explained per layer', caption: 'R² of all factors together. A layer with low total R² keeps most of its variation private or noisy.' });
    dd.draw(body => C.barsH(body, { items: r.views.map((v, i) => ({ label: v.short, value: G.r2[i].total, color: VIEW_COLORS[i] })), xLabel: 'R² (all factors)', xFormat: '.0%', domain: [0, 1], valueFormat: x => pct(x), labelWidth: 48, rowH: 16 }));
    const assoc = r.assoc.filter(x => x.key !== 'site' || S.ds.samples.length < 60);
    const e = panel(grid, { title: 'What does each factor track?', wide: false, caption: `−log₁₀ p for the association between factor scores and each covariate (ANOVA for categories, Spearman for numbers). Tumour-only covariates are tested within tumours. Cell label: η² or ρ. Covariates were never shown to the model; a strong cell is a discovery or a warning (batch).` });
    e.draw(body => C.heatmap(body, { width: 'onehalf', rows: assoc.map(x => x.label), cols: d3.range(K).map(k => `F${k + 1}`), values: assoc.map(x => x.perFactor.map(f => (isFinite(f.p) ? Math.min(30, -Math.log10(Math.max(f.p, 1e-300))) : null))), scale: 'sequential', domain: [0, 10],
      cellLabels: (v, cc) => { const f = assoc[cc.i].perFactor[cc.j]; return isFinite(f.stat) ? (f.kind === 'η²' ? f.stat.toFixed(2) : f.stat.toFixed(2)) : ''; }, cell: 34, cellH: 16, ml: 112, mb: 20, colAngle: 0, legendTitle: '−log₁₀ p',
      tooltip: cc => { const f = assoc[cc.i].perFactor[cc.j]; return `<b>${esc(assoc[cc.i].label)}</b> · Factor ${cc.j + 1}<br>${f.kind} = ${fmt(f.stat)}<br>p = ${fmtP(f.p)} · n = ${f.n}`; } }));
    sec(root, 'Inspect one factor', 'Scores show which samples sit high or low; weights show which features build the factor in each layer.');
    const ctr2 = h('div', { class: 'controls' }); ctr2.append(seg(d3.range(K).map(k => ({ value: k, label: `F${k + 1}` })), S.factor, v => { S.factor = +v; render(); })); root.append(ctr2);
    const g2 = h('div', { class: 'grid three' }); root.append(g2);
    const f = panel(g2, { title: `Factor ${S.factor + 1} scores by group`, caption: 'Each point is a sample. Box: median and IQR.' });
    f.draw(body => C.strip(body, { categories: contrast().levels.map((l, i) => ({ name: l, label: contrast().labels[i], color: GROUP_COLORS[i] })), points: r.samples.map((id, i) => ({ cat: r.groups[i], y: G.Z[i][S.factor], label: id })), yLabel: `Factor ${S.factor + 1} score` }));
    r.views.forEach((v, m) => {
      const W = G.W[m].map((w, rr) => ({ label: v.modelLabels[rr], value: w[S.factor] })).sort((x, y) => Math.abs(y.value) - Math.abs(x.value)).slice(0, 12).sort((x, y) => y.value - x.value);
      const p = panel(g2, { title: `${v.short}: top weights`, caption: `Twelve largest absolute weights (R² = ${pct(G.r2[m].perFactor[S.factor])}).` });
      p.draw(body => C.lollipop(body, { items: W, color: VIEW_COLORS[m], xLabel: 'Weight', labelWidth: 78, maxChars: 20, italic: isHuman() && m < 2 }));
    });
    sec(root, 'Do the layers agree, and how many factors?', 'Two quick checks on the model.');
    const g3 = h('div', { class: 'grid' }); root.append(g3);
    const ag = panel(g3, { title: 'Agreement between layers', caption: 'Spearman correlation between the sample-distance matrices of two layers (a Mantel statistic). High values mean the layers place samples similarly. No p-value is shown: distances are not independent.' });
    ag.draw(body => { const n = r.views.length; const M = d3.range(n).map(i => d3.range(n).map(j => (i === j ? 1 : (r.agree.find(x => (x.a === i && x.b === j) || (x.a === j && x.b === i)) || {}).r))); C.heatmap(body, { rows: r.views.map(v => v.short), cols: r.views.map(v => v.short), values: M, scale: 'sequential', domain: [0, 1], cellLabels: v => fmt(v), cell: 34, cellH: 24, ml: 56, mb: 34, colAngle: 0, legendTitle: 'Mantel ρ' }); });
    const kx = panel(g3, { title: 'More factors: what changes?', caption: 'Refit the same inputs with K = 2…8 and compare total R² per layer. The comparison does not pick an "optimal" K; extra factors capture smaller, often layer-specific structure.' });
    const run = h('button', { class: 'primary', onclick: async () => { run.disabled = true; kx.body.classList.add('loading'); try { S.kexp = await call('kexp', { Ks: [2, 3, 4, 5, 6, 7, 8] }); } finally { run.disabled = false; kx.body.classList.remove('loading'); } drawK(); } }, 'Compare K = 2 to 8');
    kx.el.insertBefore(h('div', { class: 'controls' }), kx.body).append(run);
    function drawK() { if (!S.kexp) { kx.body.innerHTML = '<div class="empty">Takes a few seconds: seven model fits.</div>'; return; }
      kx.draw(body => C.lines(body, { series: r.views.map((v, m) => ({ label: v.short, color: VIEW_COLORS[m], dots: true, points: S.kexp.map(x => [x.K, x.r2[m].total]) })), xLabel: 'Number of factors (K)', yLabel: 'Total R²', yFormat: '.0%', xTicks: 7, yDomain: [0, 1] })); }
    drawK();
    if (S.ds.batch) pause(root, S.params.batchAdjust ? 'Batch adjustment is on. Which factors appeared that were hidden before, and what do they track now?' : `Factor 1 explains most of the variance in every tissue. Look at its row in the association heatmap: what does it track? Then turn on “Batch adjustment” in the analysis settings.`,
      '<p>Without adjustment the first factor follows <b>year</b>, which here means a different season, a different mass-spectrometry plate and a separate database search. It is shared by all three tissues because the batch touched every sample of that year. Centring every protein within year removes this location shift; afterwards the factors can pick up the timepoint, colony and pathogen-related structure. Notice that the <i>paired</i> tests barely changed: both samples of a pair come from the same year, so pairing had already cancelled the batch. Batch adjustment is a modelling choice with assumptions (balanced design within batch); report it.</p>');
    pause(root, 'Which factor carries information from more than one layer? Is it the one that tracks your group?', '<p>Read the R² heatmap row-wise and the association heatmap column-wise. A factor with R² in all layers and a strong association with the group is the shared biological axis. Factors active in a single layer often track technical structure (a batch, a sample-preparation effect) or layer-specific biology (e.g. phosphorylation without protein change). Factors that track sex are a useful positive control: they should appear even though sex was never given to the model.</p>');
  };

  RENDER.network = root => {
    const K = S.res.gfa.K, q = S.netQ || (S.netQ = { factor: S.factor, mode: 'factor', per: 8, minRho: 0.6 });
    const ctr = h('div', { class: 'controls' });
    ctr.append(select('Factor', d3.range(K).map(k => ({ value: k, label: `Factor ${k + 1}` })), q.factor, v => { q.factor = +v; draw(); }),
      select('Links show', [{ value: 'factor', label: 'Factor → feature weights' }, { value: 'corr', label: 'Cross-layer correlation' }], q.mode, v => { q.mode = v; draw(); }),
      select('Features per layer', [4, 6, 8, 10, 12].map(v => ({ value: v, label: v })), q.per, v => { q.per = +v; draw(); }),
      select('Minimum |ρ|', [0.4, 0.5, 0.6, 0.7, 0.8].map(v => ({ value: v, label: v.toFixed(1) })), q.minRho, v => { q.minRho = +v; draw(); }));
    root.append(ctr);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const a = panel(grid, { title: 'Molecular network', wide: true, caption: '' });
    const b = panel(grid, { title: 'Selected feature', caption: 'Click a node to show its model-input values against the factor scores.' });
    b.body.innerHTML = '<div class="empty">Click a feature node.</div>';
    const legendItems = S.ds.views.map((v, i) => ({ label: v.short, color: VIEW_COLORS[i] })).concat([{ label: 'positive', color: '#b2182b', type: 'line' }, { label: 'negative', color: '#2166ac', type: 'line' }]);
    async function draw() {
      a.body.classList.add('loading');
      try {
        const net = await call('network', { q: Object.assign({}, q, { within: false }) });
        a.draw(body => C.network(body, { width: 'onehalf', aspect: 0.78, nodes: net.nodes.map(n => Object.assign(n, { color: n.view >= 0 ? VIEW_COLORS[n.view] : '#000' })), links: net.links, posColor: '#b2182b', negColor: '#2166ac', legendItems,
          nodeTooltip: n => n.type === 'factor' ? `Factor ${q.factor + 1}` : `<b>${esc(n.label)}</b><br>${esc(S.ds.views[n.view].name)}<br>weight ${fmt(n.w, 3)}`,
          linkTooltip: l => l.kind === 'factor' ? `weight ${fmt(l.value, 3)}` : `ρ = ${fmt(l.value)} · n = ${l.n}`,
          onClick: n => { if (n.type !== 'factor') showNode(n); } }));
        a.setCaption(q.mode === 'factor' ? `The ${q.per} largest absolute weights of Factor ${q.factor + 1} in each layer. Edge width ∝ |weight|, scaled across layers; red positive, blue negative.` : `Spearman correlations (|ρ| ≥ ${q.minRho}, ≥ 10 observed pairs) between the top-weighted features of Factor ${q.factor + 1} in different layers; up to 120 strongest edges. Co-variation is not a physical interaction and has no direction.`);
      } catch (e) { a.body.innerHTML = `<div class="empty">${esc(e.message)}</div>`; } finally { a.body.classList.remove('loading'); }
    }
    function showNode(n) {
      const v = S.ds.views[n.view], ids = S.res.samples; const vals = ids.map(id => { const x = v.values[v.sample_ids.indexOf(id)][n.j]; return x == null ? NaN : x; });
      const pts = ids.map((id, i) => ({ id, x: S.res.gfa.Z[i][q.factor], y: vals[i], group: S.res.groups[i], shape: SHAPES[S.res.codes[i]], label: id }));
      b.el.querySelector('.panel-title').textContent = `${n.label} · ${v.short}`;
      b.draw(body => C.scatter(body, { points: pts, groups: groupsSpec(), xLabel: `Factor ${q.factor + 1} score`, yLabel: v.unit.length > 24 ? 'log₂ abundance' : v.unit, fit: true, stat: `ρ = ${fmt(E.spearman(pts.map(p => p.x), pts.map(p => p.y)).r)}` }));
    }
    draw();
  };

  RENDER.pathways = root => {
    const K = S.res.gfa.K, cats = S.gs ? Object.keys(S.gs.catalogs) : [];
    const q = S.pathQ || (S.pathQ = { view: 0, source: 'factor', factor: S.factor, method: 'rank', sign: 'pos', top: 40, catalog: cats.includes('hallmark') ? 'hallmark' : cats[0] });
    if (!S.gs) { root.append(h('div', { class: 'empty' }, 'Loading gene-set catalogues…')); setTimeout(() => S.chapter === 'pathways' && render(), 800); return; }
    const ctr = h('div', { class: 'controls' });
    ctr.append(select('Layer', S.ds.views.map((v, i) => ({ value: i, label: v.name })), q.view, v => { q.view = +v; draw(); }),
      select('Score features by', [{ value: 'factor', label: 'Factor weights' }, { value: 'de', label: 'Differential t statistic' }], q.source, v => { q.source = v; render(); }),
      ...(q.source === 'factor' ? [select('Factor', d3.range(K).map(k => ({ value: k, label: `Factor ${k + 1}` })), q.factor, v => { q.factor = +v; draw(); })] : []),
      select('Catalogue', cats.map(c => ({ value: c, label: S.gs.labels[c] })), q.catalog, v => { q.catalog = v; draw(); }),
      select('Test', [{ value: 'rank', label: 'Rank-based (all scored features)' }, { value: 'ora', label: 'Over-representation (top features)' }], q.method, v => { q.method = v; render(); }),
      ...(q.method === 'ora' ? [select('Direction', [{ value: 'pos', label: 'Positive' }, { value: 'neg', label: 'Negative' }, { value: 'abs', label: 'Largest |score|' }], q.sign, v => { q.sign = v; draw(); }), select('Top features', [20, 40, 60, 100].map(v => ({ value: v, label: v })), q.top, v => { q.top = +v; draw(); })] : []));
    root.append(ctr);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const a = panel(grid, { title: 'Enriched gene sets', wide: true, caption: '' });
    const t = panel(grid, { title: 'All results', wide: true, caption: 'Sorted by p. Download the table for the full list.' });
    async function draw() {
      a.body.classList.add('loading');
      try {
        const res = await call('enrich', { q: Object.assign({}, q) }); S.enrichRes = res;
        const top = res.results.slice(0, 15);
        if (!top.length) { a.body.innerHTML = '<div class="empty">No gene set passes the size limits for this query. Try more model features, the differential statistic, or another catalogue.</div>'; t.body.innerHTML = ''; return; }
        const isRank = res.method === 'rank';
        a.draw(body => C.dotplot(body, { width: 'onehalf', items: top.map(x => ({ id: x.id, name: x.name, x: isRank ? x.auc : x.ratio, size: isRank ? x.size : x.hits, color: -Math.log10(Math.max(x.q, 1e-30)), raw: x })), xLabel: isRank ? 'AUC (set members vs all other features)' : 'Fold enrichment', vline: isRank ? 0.5 : 1, xDomain: isRank ? [0, 1] : undefined,
          sizeTitle: isRank ? 'Set size' : 'Hits', labelWidth: 150, maxChars: 40,
          tooltip: d => `<b>${esc(d.name)}</b><br>${isRank ? `AUC ${fmt(d.raw.auc)} · z ${fmt(d.raw.z)}` : `${d.raw.hits}/${d.raw.size} · ×${fmt(d.raw.ratio)}`}<br>p = ${fmtP(d.raw.p)} · q = ${fmtP(d.raw.q)}<br><i>${esc(d.raw.genes.slice(0, 6).map(g => labelOf(q.view, g)).join(', '))}</i>` }));
        a.setCaption(isRank ? `Mann–Whitney rank test of ${q.source === 'factor' ? `Factor ${q.factor + 1} weights` : 'moderated t statistics'} for each set against all other scored features in this layer (${res.nRanked} features), BH-adjusted. AUC > 0.5: members score higher (positive weights / higher in ${esc(contrast().labels[1].toLowerCase())}); &lt; 0.5: lower. Background = features actually measured, not the genome.`
          : `Hypergeometric test of the top ${q.top} ${q.sign === 'neg' ? 'negative' : q.sign === 'pos' ? 'positive' : 'absolute'} scores against the ${res.nUniverse} scored features of this layer, BH-adjusted. Fold enrichment = observed / expected hits.`);
        t.body.innerHTML = `<div class="table-wrap"><table class="data"><thead><tr><th>Gene set</th><th class="num">${isRank ? 'Size' : 'Hits/size'}</th><th class="num">${isRank ? 'AUC' : 'Fold'}</th><th class="num">p</th><th class="num">q</th><th>Leading features</th></tr></thead><tbody>${res.results.slice(0, 80).map(x => `<tr><td>${esc(x.name)}</td><td class="num">${isRank ? x.size : `${x.hits}/${x.size}`}</td><td class="num">${fmt(isRank ? x.auc : x.ratio)}</td><td class="num">${fmtP(x.p)}</td><td class="num">${fmtP(x.q)}</td><td><i>${esc(x.genes.slice(0, 6).map(g => labelOf(q.view, g)).join(', '))}</i></td></tr>`).join('')}</tbody></table></div>`;
        const tools = t.el.querySelector('.panel-tools'); tools.querySelectorAll('.csvbtn').forEach(x => x.remove());
        tools.prepend(h('button', { class: 'csvbtn', onclick: () => download(`${S.ds.id}_enrichment_${q.catalog}.csv`, csv(res.results, [['set_id', 'id'], ['name', 'name'], ['size', 'size'], ['hits', 'hits'], ['auc', 'auc'], ['fold', 'ratio'], ['p', 'p'], ['q_BH', 'q'], ['features', r => r.genes.join(';')]])) }, 'CSV'));
      } catch (e) { a.body.innerHTML = `<div class="empty">${esc(e.message)}</div>`; } finally { a.body.classList.remove('loading'); }
    }
    draw();
    pause(root, 'Run the same query as a rank test and as over-representation of the top 40. Why can the answers differ?', '<p>Over-representation depends on an arbitrary cut-off and ignores everything below it; it also asks only about one direction. The rank test uses every scored feature and is less sensitive to the cut-off, but it answers "do members tend to score higher?", which large sets can satisfy with small shifts. Both depend on the background: only features that were measured, passed filtering, and (for factor weights) entered the model.</p>');
    if (!isHuman()) root.append(h('p', { class: 'note' }, 'Honey bee GO annotations come from UniProtKB entries matched to each RefSeq protein and are mostly electronic (inferred from sequence). Treat enrichments as hypotheses.'));
  };
  function labelOf(m, sym) { const v = S.ds.views[m]; const j = v.symbols.indexOf(sym); return j >= 0 ? v.labels[j] : sym; }

  RENDER.predict = root => {
    const q = S.predQ || (S.predQ = { target: '', topN: 20 });
    const targets = [{ value: '', label: `${contrast().labels[1]} vs ${contrast().labels[0]}` }].concat(S.ds.covariates.filter(c => c.type === 'categorical' && c.key !== contrast().key && c.key !== 'site').map(c => ({ value: c.key, label: c.label })));
    const ctr = h('div', { class: 'controls' });
    ctr.append(select('Predict', targets, q.target, v => { q.target = v; S.pred = null; render(); }), select('Features selected per fold', [5, 10, 20, 50].map(v => ({ value: v, label: v })), q.topN, v => { q.topN = +v; S.pred = null; render(); }));
    const go = h('button', { class: 'primary', onclick: async () => { go.disabled = true; setStatus('busy', 'Cross-validating'); try { S.pred = await call('predict', { q: Object.assign({}, q, { target: q.target || null }) }); setStatus('ready', 'Analysis ready'); } catch (e) { setStatus('error', e.message); } go.disabled = false; render(); } }, S.pred ? 'Run again' : 'Run cross-validation');
    ctr.append(go); root.append(ctr);
    root.append(h('div', { class: 'prose' }, `<p>Each layer gets the same classifier: pick the ${q.topN} features with the largest |t| in the training samples, standardise with training statistics, fit a ridge-penalised logistic regression, and score the held-out samples. Five-fold stratified cross-validation, repeated five times. "All layers" concatenates the layers before selection (early integration). The <b>leaky</b> pipeline picks features once on all samples before splitting, which is a very common mistake.</p>`));
    if (!S.pred) { root.append(h('div', { class: 'empty' }, 'Choose a target and run the cross-validation (a few seconds).')); return; }
    const P = S.pred, names = P.honest.map(x => x.name), col = i => (i < S.ds.views.length ? VIEW_COLORS[i] : '#16191d');
    stats(root, [[P.n, 'samples with a label'], [`${P.nPos} / ${P.n - P.nPos}`, `positive / negative (${esc(P.target)})`], [fmt(d3.max(P.honest, x => x.auc)), 'best honest AUC'], [fmt(d3.max(P.permL, x => x.auc)), 'leaky AUC with shuffled labels']]);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const a = panel(grid, { title: 'Honest vs leaky cross-validation', caption: 'Mean AUC ± SD over five repeats of five-fold CV. Coloured: honest pipeline (selection inside each training fold); grey: leaky pipeline (selection on all samples first). Dashed line: chance (0.5).' });
    a.draw(body => C.barsH(body, { items: names.flatMap((nm, i) => [{ label: `${nm} · honest`, value: P.honest[i].auc, err: P.honest[i].aucSD, color: col(i) }, { label: `${nm} · leaky`, value: P.leaky[i].auc, err: P.leaky[i].aucSD, color: '#bdbdbd' }]), domain: [0, 1], vline: 0.5, xLabel: 'Cross-validated AUC', labelWidth: 92, rowH: 11 }));
    const b = panel(grid, { title: 'ROC curves (honest pipeline)', caption: 'Receiver operating characteristic from scores averaged over repeats.' });
    b.draw(body => C.lines(body, { series: P.honest.map((x, i) => ({ label: `${x.name} (${fmt(x.auc)})`, color: col(i), points: x.roc, width: i === names.length - 1 ? 1.4 : 1 })), xLabel: 'False-positive rate', yLabel: 'True-positive rate', xDomain: [0, 1], yDomain: [0, 1], diagonal: true, step: true, aspect: 0.86 }));
    const c = panel(grid, { title: 'The same pipelines on shuffled labels', wide: true, caption: 'Labels permuted once, so there is nothing to learn. The honest pipeline stays near 0.5; selecting features on all samples first finds noise that "predicts" the held-out samples, because those samples helped pick the features.' });
    c.draw(body => C.barsH(body, { width: 'onehalf', items: names.flatMap((nm, i) => [{ label: `${nm} · honest`, value: P.permH[i].auc, err: P.permH[i].aucSD, color: col(i) }, { label: `${nm} · leaky`, value: P.permL[i].auc, err: P.permL[i].aucSD, color: '#bdbdbd' }]), domain: [0, 1], vline: 0.5, xLabel: 'Cross-validated AUC (shuffled labels)', labelWidth: 92, rowH: 10 }));
    pause(root, 'Does combining the layers always predict better than the best single layer?', '<p>Not necessarily. Concatenation adds many uninformative features and the selection step must find the useful ones among more candidates; with few samples that costs as much as the extra information gains. Integration helps most when layers carry complementary signal, and its benefit must be shown with an honest, nested evaluation like this one.</p>');
  };

  RENDER.challenge = root => {
    const scen = [
      { id: 'none', title: 'Original pairing', text: 'The matched data, aligned by sample ID.' },
      { id: 'shuffle', title: 'Shuffle rows, keep their IDs', text: 'Every layer’s rows are reordered together with their IDs. Alignment by ID should recover the original result exactly.' },
      { id: 'mispair', title: `Assign the ${S.ds.views[S.ds.views.length - 1].short.toLowerCase()} layer to the wrong samples`, text: 'Values move, labels stay: a sample-tracking error.' },
      { id: 'mask', title: 'Hide 10% of measurements', text: 'Random values removed from every layer (a teaching perturbation).' },
      { id: 'subsample', title: `Keep half of the ${S.ds.pairing ? S.ds.pairing.label.toLowerCase() + 's' : 'samples'}`, text: 'Same design, half the sample size.' }];
    const cur = S.challengeScen || 'none';
    const cards = h('div', { class: 'cards' });
    scen.forEach(sc => { const c = h('div', { class: 'ds-card' + (sc.id === cur ? ' active' : ''), onclick: () => runScenario(sc.id) }, `<h3>${esc(sc.title)}</h3><p>${esc(sc.text)}</p>`); cards.append(c); });
    root.append(h('p', { class: 'prose' }, 'Predict what will happen, then choose a scenario. The source data stay unchanged; the scenario is re-run through the same pipeline with your current settings. Other chapters keep showing the original analysis.'));
    root.append(cards);
    const base = S.res, R = S.challengeRes && S.challengeScen === cur ? S.challengeRes : (cur === 'none' ? base : null);
    if (!R) { root.append(h('div', { class: 'empty' }, 'Running…')); return; }
    const agree = x => d3.mean(x.agree, a => a.r), nsig = x => x.views.map(v => v.de.filter(r => r.q < sigQ).length);
    stats(root, [[fmt(agree(R)), 'mean agreement between layers (Mantel ρ)'], [fmt(agree(base)), 'original pairing'], ...R.views.map((v, i) => [`${nsig(R)[i]}`, `${esc(v.short)} q < 0.05 (original ${nsig(base)[i]})`])]);
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const a = panel(grid, { title: 'Joint PCA in this scenario', caption: 'Same settings as the main analysis.' });
    a.draw(body => { const key = contrast().key; const pts = R.samples.map((id, i) => ({ id, x: R.joint.scores[0][i], y: R.joint.scores[1][i], group: R.groups[i], shape: SHAPES[R.codes[i]], label: id })); C.scatter(body, { points: pts, groups: groupsSpec(), xLabel: `PC1 (${pct(R.joint.varExp[0])})`, yLabel: `PC2 (${pct(R.joint.varExp[1])})`, ellipses: true, zeroLines: true }); return key; });
    const b = panel(grid, { title: 'Variance explained: scenario vs original', caption: 'Total R² of the factor model per layer.' });
    b.draw(body => C.barsH(body, { items: R.views.flatMap((v, i) => [{ label: `${v.short} · original`, value: base.gfa.r2[i].total, color: '#c4c4c4' }, { label: `${v.short} · scenario`, value: R.gfa.r2[i].total, color: VIEW_COLORS[i] }]), domain: [0, 1], xLabel: 'Total R²', xFormat: '.0%', valueFormat: x => pct(x), labelWidth: 86, rowH: 12 }));
    const c = panel(grid, { title: 'Agreement between each pair of layers', caption: 'Mantel ρ between sample-distance matrices. A mis-assigned layer loses its agreement with the others; the remaining pair is unaffected.' });
    c.draw(body => { const items = []; R.agree.forEach((x, k) => { const nm = `${R.views[x.a].short}–${R.views[x.b].short}`; items.push({ label: `${nm} · original`, value: base.agree[k].r, color: '#c4c4c4' }, { label: `${nm} · scenario`, value: x.r, color: '#16191d' }); }); C.barsH(body, { items, domain: [Math.min(0, d3.min(items, d => d.value)), 1], xLabel: 'Mantel ρ', labelWidth: 112, rowH: 11 }); });
    const st = panel(grid, { title: 'How stable are the factors?', caption: 'Refit the original model on random 80% subsets of samples (8 times) and correlate each refitted factor with the matching original factor on the shared samples. Stable factors stay near 1; unstable ones are better not interpreted.' });
    const sb = h('button', { class: 'primary', onclick: async () => { sb.disabled = true; st.body.classList.add('loading'); try { S.stab = await call('stability', { q: { B: 8, frac: 0.8 } }); } finally { sb.disabled = false; st.body.classList.remove('loading'); } drawStab(); } }, 'Run stability check');
    st.el.insertBefore(h('div', { class: 'controls' }), st.body).append(sb);
    function drawStab() { if (!S.stab) { st.body.innerHTML = '<div class="empty">Eight refits; about ten seconds.</div>'; return; }
      st.draw(body => C.strip(body, { categories: S.stab.perFactor.map((_, k) => ({ name: `F${k + 1}`, label: `F${k + 1}`, color: '#16191d' })), points: S.stab.perFactor.flatMap((v, k) => v.map(x => ({ cat: `F${k + 1}`, y: x, label: `F${k + 1}` }))), yLabel: '|r| with original factor', tooltip: p => `${p.cat}: |r| = ${fmt(p.y)}` })); }
    drawStab();
    pause(root, 'Which is more damaging: a different row order, or a wrong sample identity? Explain the difference.', '<p>Row order is harmless when every step aligns layers by sample ID, which is why "shuffle rows" reproduces the original numbers exactly. A wrong identity silently breaks the link between layers: single-layer tests are still valid, but the cross-layer agreement collapses and shared factors lose the mis-assigned layer. Nothing in the data warns you, except this kind of agreement check. Hiding values or halving the sample size mostly costs power: fewer discoveries, less stable factors.</p>');
    async function runScenario(id) {
      S.challengeScen = id; if (id === 'none') { S.challengeRes = null; render(); return; }
      S.challengeRes = null; render(); setStatus('busy', 'Running scenario');
      try { S.challengeRes = await call('run', { params: Object.assign({}, S.params, { challenge: id }) }); setStatus('ready', 'Scenario ready'); } catch (e) { setStatus('error', e.message); }
      if (S.chapter === 'challenge') render();
    }
    if (cur !== 'none' && !S.challengeRes) runScenario(cur);
  };

  RENDER.export = root => {
    const p = S.params;
    root.append(h('div', { class: 'prose' }, `<p>Every figure in this lab is drawn at journal size: <b>89 mm</b> (single column), <b>120 mm</b> or <b>183 mm</b> (double column), with 6–7 pt Helvetica/Arial text and 0.6 pt axes. Use <span class="kbd">SVG</span> on any panel for an editable vector file (Illustrator, Inkscape, Affinity) or <span class="kbd">PNG</span> for a 600 dpi raster. The bundle below collects the tables, the settings and every figure you have opened in this session.</p>`));
    const grid = h('div', { class: 'grid' }); root.append(grid);
    const a = panel(grid, { title: 'Results bundle', caption: `${S.figs.size} figure${S.figs.size === 1 ? '' : 's'} opened so far will be included as SVG.` });
    a.body.innerHTML = `<p class="note">Settings, sample sheet, differential results per layer, factor scores and weights, variance explained and covariate associations, as CSV and JSON.</p>`;
    a.body.append(h('button', { class: 'primary', onclick: bundle }, 'Download results bundle (.zip)'));
    const b = panel(grid, { title: 'Companion R walkthrough', caption: 'The same pipeline with MOFA2, limma and fgsea, using the matrices shipped with this site.' });
    b.body.innerHTML = `<p class="note">Open <a href="r/walkthrough.R">r/walkthrough.R</a> in RStudio with the folder <a href="r/data/">r/data/</a>. Every chapter also has its own snippet behind “R code for this chapter”.</p>`;
    const c = panel(grid, { title: 'Current settings', wide: true, caption: 'Recorded in the bundle as settings.json.' });
    c.body.innerHTML = `<table class="data"><tbody>${Object.entries(Object.assign({ dataset: S.ds.id }, p)).map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(JSON.stringify(v))}</td></tr>`).join('')}</tbody></table>`;
    pause(root, 'Choose one signal from each layer. What experiment would you run next?', '<p>Good next questions are specific and testable: a factor that tracks a covariate suggests a stratified re-analysis; a feature that is strong in protein but absent in RNA points to post-transcriptional regulation; a tissue-wide protein change in bees suggests a colony-level exposure that can be measured directly (residues, pathogen loads, pollen sources).</p>');
  };

  async function bundle() {
    const zip = new JSZip(), r = S.res, d = S.ds;
    zip.file('settings.json', JSON.stringify(Object.assign({ dataset: d.id, generated: new Date().toISOString() }, S.params), null, 2));
    zip.file('samples.csv', csv(d.samples, Object.keys(d.samples[0]).map(k => [k, k])));
    r.views.forEach((v, m) => {
      zip.file(`differential_${d.views[m].key}.csv`, csv(v.de, [['feature', 'f'], ['label', 'label'], ['difference', 'diff'], ['moderated_t', 't'], ['p', 'p'], ['q_BH', 'q'], ['welch_p', 'welchP'], ['welch_q', 'welchQ']]));
      zip.file(`factor_weights_${d.views[m].key}.csv`, csv(r.gfa.W[m].map((w, i) => ({ f: d.views[m].features[v.modelFeatures[i]], label: v.modelLabels[i], w })), [['feature', 'f'], ['label', 'label'], ...d3.range(r.gfa.K).map(k => [`Factor${k + 1}`, x => x.w[k]])]));
    });
    zip.file('factor_scores.csv', csv(r.samples.map((id, i) => ({ id, z: r.gfa.Z[i] })), [['sample', 'id'], ...d3.range(r.gfa.K).map(k => [`Factor${k + 1}`, x => x.z[k]])]));
    zip.file('variance_explained.csv', csv(r.views.map((v, m) => ({ v: v.name, r: r.gfa.r2[m] })), [['layer', 'v'], ['total', x => x.r.total], ...d3.range(r.gfa.K).map(k => [`Factor${k + 1}`, x => x.r.perFactor[k]])]));
    zip.file('factor_covariate_association.csv', csv(r.assoc.flatMap(a => a.perFactor.map((f, k) => ({ cov: a.label, k: k + 1, f }))), [['covariate', 'cov'], ['factor', 'k'], ['statistic', x => x.f.kind], ['value', x => x.f.stat], ['p', x => x.f.p], ['n', x => x.f.n]]));
    const figs = zip.folder('figures'); S.figs.forEach((fn, name) => { try { figs.file(name + '.svg', fn()); } catch (e) { /* figure no longer in the page */ } });
    zip.file('README.txt', `Multi-omics Lab (BIOC 557) results bundle\nDataset: ${d.title}\nSources: ${d.source.papers.map(x => x.cite + ' doi:' + x.doi).join(' | ')}\nFigures are SVG at physical size (mm).\n`);
    const blob = await zip.generateAsync({ type: 'blob' });
    const a = h('a', { href: URL.createObjectURL(blob), download: `${d.id}_results.zip` }); document.body.append(a); a.click(); setTimeout(() => a.remove(), 500);
  }

  // ---------- boot ----------
  $('settingsBtn').onclick = () => { const s = $('settings'); s.hidden = !s.hidden; $('settingsBtn').setAttribute('aria-expanded', String(!s.hidden)); };
  $('menuBtn').onclick = () => $('side').classList.toggle('open');
  $('codeClose').onclick = () => $('codeDialog').close();
  $('codeCopy').onclick = () => navigator.clipboard && navigator.clipboard.writeText($('codeBody').textContent);
  window.addEventListener('hashchange', route);
  let startKey = 'human'; try { startKey = localStorage.getItem('mol-ds') || 'human'; } catch (e) { /* ignore */ }
  const qp = new URLSearchParams(location.search).get('data'); if (qp && DATASETS.some(d => d.key === qp)) startKey = qp;
  S.chapter = (location.hash || '#study').slice(1);
  loadDataset(startKey).then(route);
})();
