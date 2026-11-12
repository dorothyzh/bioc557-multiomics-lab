/* Multi-omics teaching lab — analysis engine.
 * Pure JavaScript, no dependencies. Runs in a Web Worker (importScripts) and in Node (require).
 * Matrices are row-major Float64Array, samples x features, NaN = missing.
 */
(function (root) {
  'use strict';

  // ---------- random ----------
  function rng(seed) {
    let a = (seed >>> 0) || 1;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function shuffleIdx(n, rand) {
    const p = Array.from({ length: n }, (_, i) => i);
    for (let i = n - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); const t = p[i]; p[i] = p[j]; p[j] = t; }
    return p;
  }

  // ---------- special functions ----------
  const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61503916999185, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  function lgamma(x) {
    if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
    x -= 1; let a = LANCZOS[0]; const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += LANCZOS[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }
  function digamma(x) {
    let r = 0;
    while (x < 6) { r -= 1 / x; x += 1; }
    const f = 1 / (x * x);
    return r + Math.log(x) - 0.5 / x - f * (1 / 12 - f * (1 / 120 - f * (1 / 252 - f * (1 / 240 - f / 132))));
  }
  function trigamma(x) {
    let r = 0;
    while (x < 6) { r += 1 / (x * x); x += 1; }
    const f = 1 / (x * x);
    return r + 1 / x + f / 2 + f / x * (1 / 6 - f * (1 / 30 - f * (1 / 42 - f / 30)));
  }
  // limma::trigammaInverse
  function trigammaInverse(x) {
    if (!(x > 0)) return NaN;
    if (x > 1e7) return 1 / Math.sqrt(x);
    if (x < 1e-6) return 1 / x;
    let y = 0.5 + 1 / x;
    for (let i = 0; i < 50; i++) {
      const tri = trigamma(y);
      const dif = tri * (1 - tri / x) / psi2(y);
      y += dif;
      if (-dif / y < 1e-8) break;
    }
    return y;
  }
  function psi2(x) { // tetragamma (second derivative of digamma)
    let r = 0;
    while (x < 6) { r -= 2 / (x * x * x); x += 1; }
    const f = 1 / (x * x);
    return r - 1 / (x * x) - 1 / (x * x * x) - f * f * (0.5 - f * (1 / 6 - f * (1 / 6 - f * 0.3)));
  }
  function betacf(a, b, x) {
    const MAXIT = 300, EPS = 3e-15, FPMIN = 1e-300;
    let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
      const del = d * c; h *= del; if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  function ibeta(x, a, b) { // regularized incomplete beta I_x(a,b)
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    if (x < (a + 1) / (a + b + 2)) return bt * betacf(a, b, x) / a;
    return 1 - bt * betacf(b, a, 1 - x) / b;
  }
  // two-sided p for t with df (df may be large / Infinity)
  function pt2(t, df) {
    if (!isFinite(t)) return isNaN(t) ? NaN : 0;
    if (!isFinite(df) || df > 1e6) return 2 * pnormUpper(Math.abs(t));
    return ibeta(df / (df + t * t), df / 2, 0.5);
  }
  function pnormUpper(z) { return z >= 0 ? erfc(z / Math.SQRT2) / 2 : 1 - erfc(-z / Math.SQRT2) / 2; }
  function erfc(x) {
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 +
      t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }
  // upper tail of F distribution
  function pf(F, d1, d2) { if (!(F > 0)) return 1; return ibeta(d2 / (d2 + d1 * F), d2 / 2, d1 / 2); }
  function lchoose(n, k) { return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1); }
  // P(X >= k), X ~ Hypergeom(N population, K successes, n draws)
  function phyperUpper(k, K, N, n) {
    const hi = Math.min(K, n); if (k > hi) return 0;
    const lo = Math.max(0, n - (N - K)); k = Math.max(k, lo);
    const denom = lchoose(N, n); let s = 0;
    for (let i = k; i <= hi; i++) s += Math.exp(lchoose(K, i) + lchoose(N - K, n - i) - denom);
    return Math.min(1, s);
  }
  function bh(p) {
    const n = p.length, idx = [];
    for (let i = 0; i < n; i++) if (isFinite(p[i])) idx.push(i);
    idx.sort((a, b) => p[b] - p[a]);
    const q = new Array(n).fill(NaN); const m = idx.length; let cm = 1;
    idx.forEach((i, r) => { const rank = m - r; cm = Math.min(cm, p[i] * m / rank); q[i] = Math.min(1, cm); });
    return q;
  }

  // ---------- basic helpers ----------
  const mean = a => { let s = 0, c = 0; for (const v of a) if (v === v) { s += v; c++; } return c ? s / c : NaN; };
  function quantile(sorted, q) {
    if (!sorted.length) return NaN; const pos = (sorted.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
    return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
  }
  function median(a) { const s = a.filter(v => v === v).sort((x, y) => x - y); return quantile(s, 0.5); }
  function ranks(a) { // average ranks for ties, input without NaN
    const idx = a.map((v, i) => i).sort((i, j) => a[i] - a[j]); const r = new Array(a.length);
    for (let i = 0; i < idx.length;) { let j = i; while (j + 1 < idx.length && a[idx[j + 1]] === a[idx[i]]) j++; const avg = (i + j) / 2 + 1; for (let k = i; k <= j; k++) r[idx[k]] = avg; i = j + 1; }
    return r;
  }
  function pearson(x, y) {
    let n = 0, sx = 0, sy = 0; for (let i = 0; i < x.length; i++) if (x[i] === x[i] && y[i] === y[i]) { n++; sx += x[i]; sy += y[i]; }
    if (n < 3) return { r: NaN, n }; const mx = sx / n, my = sy / n; let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < x.length; i++) if (x[i] === x[i] && y[i] === y[i]) { const a = x[i] - mx, b = y[i] - my; sxy += a * b; sxx += a * a; syy += b * b; }
    return { r: sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : NaN, n };
  }
  function spearman(x, y) {
    const xs = [], ys = []; for (let i = 0; i < x.length; i++) if (x[i] === x[i] && y[i] === y[i]) { xs.push(x[i]); ys.push(y[i]); }
    if (xs.length < 3) return { r: NaN, n: xs.length };
    const r = pearson(ranks(xs), ranks(ys)).r; return { r, n: xs.length, p: corP(r, xs.length) };
  }
  function corP(r, n) { if (!(n > 2) || !isFinite(r)) return NaN; const t = r * Math.sqrt((n - 2) / Math.max(1e-15, 1 - r * r)); return pt2(t, n - 2); }

  // ---------- small dense linear algebra ----------
  function cholInv(A, k) { // A: Float64Array k*k SPD -> inverse, logdet
    const L = new Float64Array(k * k);
    for (let i = 0; i < k; i++) for (let j = 0; j <= i; j++) {
      let s = A[i * k + j]; for (let p = 0; p < j; p++) s -= L[i * k + p] * L[j * k + p];
      if (i === j) L[i * k + i] = Math.sqrt(Math.max(s, 1e-12)); else L[i * k + j] = s / L[j * k + j];
    }
    const Li = new Float64Array(k * k);
    for (let i = 0; i < k; i++) { Li[i * k + i] = 1 / L[i * k + i]; for (let j = 0; j < i; j++) { let s = 0; for (let p = j; p < i; p++) s -= L[i * k + p] * Li[p * k + j]; Li[i * k + j] = s / L[i * k + i]; } }
    const inv = new Float64Array(k * k);
    for (let i = 0; i < k; i++) for (let j = 0; j <= i; j++) { let s = 0; for (let p = i; p < k; p++) s += Li[p * k + i] * Li[p * k + j]; inv[i * k + j] = inv[j * k + i] = s; }
    return inv;
  }
  function solveSPD(A, b, k) { const inv = cholInv(A, k); const x = new Float64Array(k); for (let i = 0; i < k; i++) { let s = 0; for (let j = 0; j < k; j++) s += inv[i * k + j] * b[j]; x[i] = s; } return x; }
  // symmetric eigen (Jacobi). returns {values desc, vectors (column j = vector j) row-major n*n}
  function eigSym(Ain, n) {
    const A = Float64Array.from(Ain), V = new Float64Array(n * n); for (let i = 0; i < n; i++) V[i * n + i] = 1;
    for (let sweep = 0; sweep < 60; sweep++) {
      let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
      let dg = 0; for (let i = 0; i < n; i++) dg += A[i * n + i] * A[i * n + i];
      if (off < 1e-22 * Math.max(dg, 1e-300) || off < 1e-300) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        const apq = A[p * n + q]; if (Math.abs(apq) < 1e-300) continue;
        const app = A[p * n + p], aqq = A[q * n + q], theta = (aqq - app) / (2 * apq);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const akp = A[k * n + p], akq = A[k * n + q]; A[k * n + p] = c * akp - s * akq; A[k * n + q] = s * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = A[p * n + k], aqk = A[q * n + k]; A[p * n + k] = c * apk - s * aqk; A[q * n + k] = s * apk + c * aqk; }
        for (let k = 0; k < n; k++) { const vkp = V[k * n + p], vkq = V[k * n + q]; V[k * n + p] = c * vkp - s * vkq; V[k * n + q] = s * vkp + c * vkq; }
      }
    }
    const ord = Array.from({ length: n }, (_, i) => i).sort((i, j) => A[j * n + j] - A[i * n + i]);
    const values = ord.map(i => A[i * n + i]); const vectors = new Float64Array(n * n);
    ord.forEach((src, dst) => { for (let k = 0; k < n; k++) vectors[k * n + dst] = V[k * n + src]; });
    return { values, vectors };
  }

  // top-c eigenpairs of a symmetric PSD matrix by subspace iteration + Rayleigh-Ritz (vectors returned n x n layout, first c columns valid)
  function topEig(G, n, c) {
    if (n <= 60) return eigSym(G, n);
    const b = Math.min(n, c + 6), rand = rng(11); let Q = new Float64Array(n * b);
    for (let i = 0; i < n * b; i++) Q[i] = rand() - 0.5;
    const orth = (Y) => { // modified Gram-Schmidt on columns
      for (let a = 0; a < b; a++) {
        for (let p = 0; p < a; p++) { let d = 0; for (let i = 0; i < n; i++) d += Y[i * b + a] * Y[i * b + p]; for (let i = 0; i < n; i++) Y[i * b + a] -= d * Y[i * b + p]; }
        let nn = 0; for (let i = 0; i < n; i++) nn += Y[i * b + a] ** 2; nn = Math.sqrt(nn) || 1; for (let i = 0; i < n; i++) Y[i * b + a] /= nn;
      }
      return Y;
    };
    Q = orth(Q); let prev = null;
    for (let it = 0; it < 300; it++) {
      const Y = new Float64Array(n * b);
      for (let i = 0; i < n; i++) { const gi = i * n; for (let l = 0; l < n; l++) { const g = G[gi + l]; if (g === 0) continue; const ql = l * b; for (let a = 0; a < b; a++) Y[i * b + a] += g * Q[ql + a]; } }
      // Rayleigh quotients for convergence check
      const rq = []; for (let a = 0; a < c; a++) { let s = 0; for (let i = 0; i < n; i++) s += Q[i * b + a] * Y[i * b + a]; rq.push(s); }
      Q = orth(Y);
      if (prev && rq.every((v, a) => Math.abs(v - prev[a]) <= 1e-10 * Math.abs(v) + 1e-14)) break;
      prev = rq;
    }
    // Rayleigh-Ritz: T = Q^T G Q (b x b), eigen-decompose, rotate
    const GQ = new Float64Array(n * b);
    for (let i = 0; i < n; i++) for (let l = 0; l < n; l++) { const g = G[i * n + l]; if (g === 0) continue; for (let a = 0; a < b; a++) GQ[i * b + a] += g * Q[l * b + a]; }
    const T = new Float64Array(b * b); for (let a = 0; a < b; a++) for (let d = 0; d < b; d++) { let s = 0; for (let i = 0; i < n; i++) s += Q[i * b + a] * GQ[i * b + d]; T[a * b + d] = s; }
    for (let a = 0; a < b; a++) for (let d = a + 1; d < b; d++) T[a * b + d] = T[d * b + a] = (T[a * b + d] + T[d * b + a]) / 2;
    const e = eigSym(T, b); const vectors = new Float64Array(n * n);
    for (let i = 0; i < n; i++) for (let p = 0; p < c; p++) { let s = 0; for (let a = 0; a < b; a++) s += Q[i * b + a] * e.vectors[a * b + p]; vectors[i * n + p] = s; }
    return { values: e.values.slice(0, c), vectors };
  }

  // ---------- matrix helpers ----------
  function toMatrix(values) { // array of rows (null for missing) -> {X, n, d}
    const n = values.length, d = n ? values[0].length : 0, X = new Float64Array(n * d);
    for (let i = 0; i < n; i++) { const r = values[i]; for (let j = 0; j < d; j++) { const v = r[j]; X[i * d + j] = v == null ? NaN : v; } }
    return { X, n, d };
  }
  function column(X, n, d, j) { const c = new Float64Array(n); for (let i = 0; i < n; i++) c[i] = X[i * d + j]; return c; }
  function subsetCols(X, n, d, cols) { const k = cols.length, Y = new Float64Array(n * k); for (let i = 0; i < n; i++) for (let c = 0; c < k; c++) Y[i * k + c] = X[i * d + cols[c]]; return Y; }
  function subsetRows(X, n, d, rows) { const Y = new Float64Array(rows.length * d); rows.forEach((r, i) => Y.set(X.subarray(r * d, r * d + d), i * d)); return Y; }
  function colStats(X, n, d) {
    const m = new Float64Array(d), v = new Float64Array(d), miss = new Float64Array(d);
    for (let j = 0; j < d; j++) {
      let s = 0, c = 0; for (let i = 0; i < n; i++) { const x = X[i * d + j]; if (x === x) { s += x; c++; } }
      const mu = c ? s / c : NaN; let ss = 0; for (let i = 0; i < n; i++) { const x = X[i * d + j]; if (x === x) ss += (x - mu) * (x - mu); }
      m[j] = mu; v[j] = c > 1 ? ss / (c - 1) : NaN; miss[j] = 1 - c / n;
    }
    return { mean: m, variance: v, missing: miss };
  }

  // ---------- preparation ----------
  /* params: { sampleNorm: 'none'|'median', maxMissing: 0.3, modelFeatures: 200, scaling: 'unit'|'pareto'|'none',
               center: true, impute: 'median'|'min' } */
  function prepareView(view, rowsIdx, params) {
    const full = toMatrix(view.values); const n = rowsIdx.length, d0 = full.d;
    let X = subsetRows(full.X, full.n, d0, rowsIdx);
    // QC before filtering
    const sampleMissing = new Array(n).fill(0);
    for (let i = 0; i < n; i++) { let c = 0; for (let j = 0; j < d0; j++) if (X[i * d0 + j] !== X[i * d0 + j]) c++; sampleMissing[i] = c / d0; }
    // sample normalization
    let sampleShift = new Array(n).fill(0);
    if (params.sampleNorm === 'median') {
      const meds = []; for (let i = 0; i < n; i++) meds.push(median(Array.from(X.subarray(i * d0, i * d0 + d0))));
      const grand = median(meds); sampleShift = meds.map(m => m - grand);
      for (let i = 0; i < n; i++) for (let j = 0; j < d0; j++) X[i * d0 + j] -= sampleShift[i];
    }
    // optional batch adjustment: centre each feature within batch, then restore the overall mean
    if (params.batchCodes) {
      const lv = [...new Set(params.batchCodes)];
      for (let j = 0; j < d0; j++) {
        let gs = 0, gc = 0; const bs = new Map();
        for (let i = 0; i < n; i++) { const x = X[i * d0 + j]; if (x !== x) continue; gs += x; gc++; const b = params.batchCodes[i]; const e = bs.get(b) || [0, 0]; e[0] += x; e[1]++; bs.set(b, e); }
        if (!gc) continue; const gm = gs / gc;
        for (let i = 0; i < n; i++) { const x = X[i * d0 + j]; if (x !== x) continue; const e = bs.get(params.batchCodes[i]); X[i * d0 + j] = x - e[0] / e[1] + gm; }
      }
    }
    const st0 = colStats(X, n, d0);
    // missingness + variance filter (label-free)
    const keep = []; for (let j = 0; j < d0; j++) if (st0.missing[j] <= params.maxMissing + 1e-9 && st0.variance[j] > 1e-10) keep.push(j);
    const byVar = keep.slice().sort((a, b) => st0.variance[b] - st0.variance[a]);
    const modelCols = byVar.slice(0, Math.min(params.modelFeatures, byVar.length));
    // sample distribution summaries after normalization
    const sampleBox = []; for (let i = 0; i < n; i++) {
      const s = Array.from(X.subarray(i * d0, i * d0 + d0)).filter(v => v === v).sort((a, b) => a - b);
      sampleBox.push([0.05, 0.25, 0.5, 0.75, 0.95].map(q => quantile(s, q)));
    }
    // model matrix: center + scale, keep NaN
    const k = modelCols.length, M = subsetCols(X, n, d0, modelCols), st = colStats(M, n, k);
    const scale = new Float64Array(k);
    for (let j = 0; j < k; j++) {
      const sd = Math.sqrt(st.variance[j]) || 1;
      scale[j] = params.scaling === 'unit' ? sd : params.scaling === 'pareto' ? Math.sqrt(sd) : 1;
      for (let i = 0; i < n; i++) { const x = M[i * k + j]; if (x === x) M[i * k + j] = ((params.center !== false) ? x - st.mean[j] : x) / scale[j]; }
    }
    return {
      name: view.name, unit: view.unit, n, dInput: d0, X, d: d0, keep, modelCols, M, k,
      featureMean: st0.mean, featureMissing: st0.missing, featureVar: st0.variance,
      sampleMissing, sampleBox, sampleShift,
      missingTotal: mean(Array.from(st0.missing)),
      nMissingInModel: countNaN(M)
    };
  }
  function countNaN(M) { let c = 0; for (let i = 0; i < M.length; i++) if (M[i] !== M[i]) c++; return c; }
  function imputeColumns(M, n, k, mode) { // returns copy, NaN -> column median (0 after centering) or half-min
    const Y = Float64Array.from(M);
    for (let j = 0; j < k; j++) {
      const col = []; for (let i = 0; i < n; i++) { const v = Y[i * k + j]; if (v === v) col.push(v); }
      if (!col.length) { for (let i = 0; i < n; i++) Y[i * k + j] = 0; continue; }
      col.sort((a, b) => a - b);
      const fill = mode === 'min' ? col[0] : quantile(col, 0.5);
      for (let i = 0; i < n; i++) if (Y[i * k + j] !== Y[i * k + j]) Y[i * k + j] = fill;
    }
    return Y;
  }

  // ---------- PCA ----------
  function pca(Mimp, n, k, ncomp) { // Mimp: complete, (centered) n x k
    const A = new Float64Array(n * k), mu = new Float64Array(k);
    for (let j = 0; j < k; j++) { let s = 0; for (let i = 0; i < n; i++) s += Mimp[i * k + j]; mu[j] = s / n; }
    for (let i = 0; i < n; i++) for (let j = 0; j < k; j++) A[i * k + j] = Mimp[i * k + j] - mu[j];
    const G = new Float64Array(n * n); let total = 0;
    for (let i = 0; i < n; i++) for (let l = i; l < n; l++) { let s = 0; const ri = i * k, rl = l * k; for (let j = 0; j < k; j++) s += A[ri + j] * A[rl + j]; G[i * n + l] = G[l * n + i] = s; }
    for (let i = 0; i < n; i++) total += G[i * n + i];
    const c = Math.min(ncomp, n - 1); const { values, vectors } = topEig(G, n, c);
    const scores = [], varExp = [];
    for (let p = 0; p < c; p++) { const sv = Math.sqrt(Math.max(values[p], 0)); const col = []; for (let i = 0; i < n; i++) col.push(vectors[i * n + p] * sv); scores.push(col); varExp.push(values[p] / total); }
    // loadings for top components: v = A^T u / s
    const loadings = []; for (let p = 0; p < Math.min(c, 3); p++) { const sv = Math.sqrt(Math.max(values[p], 1e-12)); const l = new Float64Array(k); for (let j = 0; j < k; j++) { let s = 0; for (let i = 0; i < n; i++) s += A[i * k + j] * vectors[i * n + p]; l[j] = s / sv; } loadings.push(Array.from(l)); }
    // deterministic sign: largest |loading| positive
    for (let p = 0; p < loadings.length; p++) { let mj = 0; loadings[p].forEach((v, j) => { if (Math.abs(v) > Math.abs(loadings[p][mj])) mj = j; }); if (loadings[p][mj] < 0) { loadings[p] = loadings[p].map(v => -v); scores[p] = scores[p].map(v => -v); } }
    return { scores, varExp, loadings };
  }

  // ---------- linear models with empirical Bayes (limma-style) ----------
  function fitFDist(s2, df) { // limma::fitFDist without covariate
    const ok = []; for (let i = 0; i < s2.length; i++) if (s2[i] > 0 && df[i] > 0 && isFinite(s2[i])) ok.push(i);
    const m = ok.length; if (m < 3) return { s02: NaN, d0: 0 };
    const e = ok.map(i => Math.log(s2[i]) - digamma(df[i] / 2) + Math.log(df[i] / 2));
    const emean = e.reduce((a, b) => a + b, 0) / m;
    let evar = e.reduce((a, b) => a + (b - emean) * (b - emean), 0) / (m - 1);
    evar -= ok.reduce((a, i) => a + trigamma(df[i] / 2), 0) / m;
    if (evar > 0) { const d0 = 2 * trigammaInverse(evar); return { d0, s02: Math.exp(emean + digamma(d0 / 2) - Math.log(d0 / 2)) }; }
    return { d0: Infinity, s02: Math.exp(emean) };
  }
  /* Per-feature OLS y ~ group (+ covariates), contrast = group[1] - group[0].
     groupCodes: 0/1 per sample. covars: array of arrays (numeric or categorical codes -> dummy-encoded by caller) */
  function featureTests(X, n, d, cols, groupCodes, covarDesign, opts) {
    opts = opts || {}; const minPerGroup = opts.minPerGroup || 3;
    const res = [];
    const p = 2 + (covarDesign ? covarDesign[0].length : 0);
    for (const j of cols) {
      const rows = []; let n0 = 0, n1 = 0;
      for (let i = 0; i < n; i++) if (X[i * d + j] === X[i * d + j]) { rows.push(i); if (groupCodes[i]) n1++; else n0++; }
      const r = { j, n0, n1, diff: NaN, t: NaN, p: NaN, s2: NaN, df: 0, welchP: NaN, welchT: NaN, mean0: NaN, mean1: NaN };
      if (n0 < minPerGroup || n1 < minPerGroup) { res.push(r); continue; }
      // group means & Welch
      let s0 = 0, s1 = 0; rows.forEach(i => { if (groupCodes[i]) s1 += X[i * d + j]; else s0 += X[i * d + j]; });
      const m0 = s0 / n0, m1 = s1 / n1; let v0 = 0, v1 = 0;
      rows.forEach(i => { const x = X[i * d + j]; if (groupCodes[i]) v1 += (x - m1) ** 2; else v0 += (x - m0) ** 2; });
      v0 /= (n0 - 1); v1 /= (n1 - 1);
      const se2 = v0 / n0 + v1 / n1, wt = (m1 - m0) / Math.sqrt(se2), wdf = se2 * se2 / ((v0 / n0) ** 2 / (n0 - 1) + (v1 / n1) ** 2 / (n1 - 1));
      r.mean0 = m0; r.mean1 = m1; r.welchT = wt; r.welchP = pt2(wt, wdf);
      if (opts.pairs) { // paired design: y ~ group + pair  <=>  moderated test on within-pair differences
        const a = new Map(), b = new Map();
        for (const i of rows) (groupCodes[i] ? b : a).set(opts.pairs[i], X[i * d + j]);
        const dif = []; a.forEach((v, k) => { if (b.has(k)) dif.push(b.get(k) - v); });
        r.nPairs = dif.length; if (dif.length < minPerGroup) { res.push(r); continue; }
        const md = dif.reduce((s, v) => s + v, 0) / dif.length, vd = dif.reduce((s, v) => s + (v - md) ** 2, 0) / (dif.length - 1);
        r.diff = md; r.s2 = vd; r.df = dif.length - 1; r.unscaled = 1 / dif.length; res.push(r); continue;
      }
      // OLS with design [1, g, covars]
      const XtX = new Float64Array(p * p), Xty = new Float64Array(p), rowv = new Float64Array(p);
      for (const i of rows) {
        rowv[0] = 1; rowv[1] = groupCodes[i]; if (covarDesign) for (let c = 0; c < covarDesign[i].length; c++) rowv[2 + c] = covarDesign[i][c];
        const y = X[i * d + j];
        for (let a = 0; a < p; a++) { Xty[a] += rowv[a] * y; for (let b = 0; b <= a; b++) XtX[a * p + b] += rowv[a] * rowv[b]; }
      }
      for (let a = 0; a < p; a++) for (let b = a + 1; b < p; b++) XtX[a * p + b] = XtX[b * p + a];
      // ridge-tiny for rank deficiency of covariates within observed rows
      for (let a = 2; a < p; a++) XtX[a * p + a] += 1e-8;
      const inv = cholInv(XtX, p); const beta = new Float64Array(p);
      for (let a = 0; a < p; a++) { let s = 0; for (let b = 0; b < p; b++) s += inv[a * p + b] * Xty[b]; beta[a] = s; }
      let rss = 0; for (const i of rows) {
        let fit = beta[0] + beta[1] * groupCodes[i]; if (covarDesign) for (let c = 0; c < covarDesign[i].length; c++) fit += beta[2 + c] * covarDesign[i][c];
        rss += (X[i * d + j] - fit) ** 2;
      }
      const df = rows.length - p; if (df < 1) { res.push(r); continue; }
      r.diff = beta[1]; r.s2 = rss / df; r.df = df; r.unscaled = inv[1 * p + 1];
      res.push(r);
    }
    // empirical Bayes moderation
    const ok = res.filter(r => r.df > 0);
    const { d0, s02 } = fitFDist(ok.map(r => r.s2), ok.map(r => r.df));
    for (const r of ok) {
      const post = isFinite(d0) ? (d0 * s02 + r.df * r.s2) / (d0 + r.df) : s02;
      const dft = isFinite(d0) ? r.df + d0 : Infinity;
      r.t = r.diff / Math.sqrt(post * r.unscaled); r.p = pt2(r.t, dft); r.dfTotal = dft; r.s2post = post;
    }
    const q = bh(res.map(r => r.p)), wq = bh(res.map(r => r.welchP));
    res.forEach((r, i) => { r.q = q[i]; r.welchQ = wq[i]; });
    return { results: res, d0, s02 };
  }

  // ---------- variance partition (type-II sum of squares per categorical covariate) ----------
  function designFromFactors(factors, n, drop) { // factors: [{key, codes[]}], returns matrix rows incl intercept
    const cols = [];
    factors.forEach((f, fi) => {
      if (fi === drop) return;
      if (f.numeric) { const m = mean(f.codes); cols.push(f.codes.map(v => (v === v ? v - m : 0))); return; }
      const lv = [...new Set(f.codes)].sort(); lv.slice(1).forEach(l => cols.push(f.codes.map(c => (c === l ? 1 : 0))));
    });
    return Array.from({ length: n }, (_, i) => [1, ...cols.map(c => c[i])]);
  }
  function rssFit(y, D) { // y with NaN, D rows -> rss over observed, tss
    const p = D[0].length, XtX = new Float64Array(p * p), Xty = new Float64Array(p); let cnt = 0, sy = 0, syy = 0;
    for (let i = 0; i < y.length; i++) { if (y[i] !== y[i]) continue; cnt++; sy += y[i]; syy += y[i] * y[i]; const r = D[i]; for (let a = 0; a < p; a++) { Xty[a] += r[a] * y[i]; for (let b = 0; b <= a; b++) XtX[a * p + b] += r[a] * r[b]; } }
    for (let a = 0; a < p; a++) { for (let b = a + 1; b < p; b++) XtX[a * p + b] = XtX[b * p + a]; if (a > 0) XtX[a * p + a] += 1e-8; }
    const beta = solveSPD(XtX, Xty, p); let rss = 0;
    for (let i = 0; i < y.length; i++) { if (y[i] !== y[i]) continue; let f = 0; for (let a = 0; a < p; a++) f += beta[a] * D[i][a]; rss += (y[i] - f) ** 2; }
    return { rss, tss: syy - sy * sy / cnt, n: cnt };
  }
  function variancePartition(X, n, d, cols, factors) {
    const full = designFromFactors(factors, n, -1), reduced = factors.map((f, fi) => designFromFactors(factors, n, fi));
    const out = factors.map(() => []), resid = [];
    for (const j of cols) {
      const y = column(X, n, d, j); const F = rssFit(y, full); if (!(F.tss > 0)) continue;
      factors.forEach((f, fi) => { const R = rssFit(y, reduced[fi]); out[fi].push(Math.max(0, (R.rss - F.rss) / F.tss)); });
      resid.push(F.rss / F.tss);
    }
    return { factors: factors.map((f, fi) => ({ key: f.key, label: f.label, values: out[fi] })), residual: resid };
  }

  // ---------- Bayesian group factor analysis (MOFA-type Gaussian model, ARD on weights) ----------
  /* views: [{M (n x k, centered/scaled, NaN missing), k}] ; returns Z, W per view, R2 */
  function gfa(views, n, K, opts) {
    opts = opts || {}; const maxIter = opts.maxIter || 400, tol = opts.tol || 3e-4, seed = opts.seed || 557;
    const a0 = 1e-3, b0 = 1e-3, KK = K * K;
    // init Z with PCA on concatenated, block-weighted, imputed data
    const blocks = views.map(v => { const Y = imputeColumns(v.M, n, v.k, 'median'); const w = 1 / Math.sqrt(v.k); for (let i = 0; i < Y.length; i++) Y[i] *= w; return { Y, k: v.k }; });
    const ktot = blocks.reduce((a, b) => a + b.k, 0), C = new Float64Array(n * ktot);
    let off = 0; for (const b of blocks) { for (let i = 0; i < n; i++) for (let j = 0; j < b.k; j++) C[i * ktot + off + j] = b.Y[i * b.k + j]; off += b.k; }
    const init = pca(C, n, ktot, K);
    const rand = rng(seed);
    const Z = new Float64Array(n * K);
    for (let k = 0; k < K; k++) { const sc = init.scores[k] || []; let ss = 0; for (let i = 0; i < n; i++) ss += (sc[i] || 0) ** 2; const sd = Math.sqrt(ss / n) || 1; for (let i = 0; i < n; i++) Z[i * K + k] = sc[i] != null ? sc[i] / sd : (rand() - 0.5) * 0.1; }
    const ZZ = new Float64Array(n * KK); // E[z z^T] per sample
    for (let i = 0; i < n; i++) for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) ZZ[i * KK + a * K + b] = Z[i * K + a] * Z[i * K + b] + (a === b ? 1e-3 : 0);
    const V = views.map(v => {
      const obs = new Uint8Array(n * v.k); let nobs = 0; for (let i = 0; i < n * v.k; i++) if (v.M[i] === v.M[i]) { obs[i] = 1; nobs++; }
      let ss = 0, c = 0; for (let i = 0; i < n * v.k; i++) if (obs[i]) { ss += v.M[i] * v.M[i]; c++; }
      const varAll = ss / Math.max(1, c);
      return { M: v.M, k: v.k, obs, nobs, W: new Float64Array(v.k * K), WW: new Float64Array(v.k * KK), alpha: new Float64Array(K).fill(1), tau: new Float64Array(v.k).fill(1 / Math.max(varAll, 1e-6)) };
    });
    // initial W by ridge regression on the PCA-initialised Z
    for (const v of V) for (let j = 0; j < v.k; j++) {
      const P = new Float64Array(KK), r = new Float64Array(K);
      for (let i = 0; i < n; i++) { const idx = i * v.k + j; if (!v.obs[idx]) continue; const y = v.M[idx]; for (let a = 0; a < K; a++) { r[a] += y * Z[i * K + a]; for (let b = 0; b < K; b++) P[a * K + b] += ZZ[i * KK + a * K + b]; } }
      for (let a = 0; a < K; a++) P[a * K + a] += 1;
      const S = cholInv(P, K);
      for (let a = 0; a < K; a++) { let s = 0; for (let b = 0; b < K; b++) s += S[a * K + b] * r[b]; v.W[j * K + a] = s; }
      for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) v.WW[j * KK + a * K + b] = v.W[j * K + a] * v.W[j * K + b] + (a === b ? 1e-3 : 0);
    }
    let iter = 0, converged = false; const elboTrace = [];
    const Zprev = new Float64Array(n * K);
    for (iter = 0; iter < maxIter; iter++) {
      Zprev.set(Z);
      // --- Z update per sample (uses current q(W), tau)
      const Ssum = new Float64Array(KK);
      for (const v of V) for (let j = 0; j < v.k; j++) { const tj = v.tau[j]; for (let t = 0; t < KK; t++) Ssum[t] += tj * v.WW[j * KK + t]; }
      for (let i = 0; i < n; i++) {
        const P = Float64Array.from(Ssum), r = new Float64Array(K);
        for (const v of V) for (let j = 0; j < v.k; j++) {
          const idx = i * v.k + j, tj = v.tau[j];
          if (!v.obs[idx]) { for (let t = 0; t < KK; t++) P[t] -= tj * v.WW[j * KK + t]; }
          else { const y = v.M[idx] * tj; for (let a = 0; a < K; a++) r[a] += y * v.W[j * K + a]; }
        }
        for (let a = 0; a < K; a++) P[a * K + a] += 1;
        const S = cholInv(P, K);
        for (let a = 0; a < K; a++) { let s = 0; for (let b = 0; b < K; b++) s += S[a * K + b] * r[b]; Z[i * K + a] = s; }
        for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) ZZ[i * KK + a * K + b] = S[a * K + b] + Z[i * K + a] * Z[i * K + b];
      }
      // --- W, alpha, tau per view/feature (sufficient statistics shared)
      const ZZsum = new Float64Array(KK); for (let i = 0; i < n; i++) for (let t = 0; t < KK; t++) ZZsum[t] += ZZ[i * KK + t];
      let sse = 0, sobs = 0;
      for (const v of V) {
        const Q = new Float64Array(KK), r = new Float64Array(K), P = new Float64Array(KK);
        for (let j = 0; j < v.k; j++) {
          Q.set(ZZsum); r.fill(0); let yy = 0, c = 0;
          for (let i = 0; i < n; i++) {
            const idx = i * v.k + j;
            if (!v.obs[idx]) { for (let t = 0; t < KK; t++) Q[t] -= ZZ[i * KK + t]; }
            else { const y = v.M[idx]; yy += y * y; c++; for (let a = 0; a < K; a++) r[a] += y * Z[i * K + a]; }
          }
          const tj = v.tau[j];
          for (let t = 0; t < KK; t++) P[t] = Q[t] * tj; for (let a = 0; a < K; a++) P[a * K + a] += v.alpha[a];
          const S = cholInv(P, K);
          for (let a = 0; a < K; a++) { let s = 0; for (let b = 0; b < K; b++) s += S[a * K + b] * r[b]; v.W[j * K + a] = s * tj; }
          for (let a = 0; a < K; a++) for (let b = 0; b < K; b++) v.WW[j * KK + a * K + b] = S[a * K + b] + v.W[j * K + a] * v.W[j * K + b];
          // expected squared error: sum y^2 - 2 w^T r + tr(E[ww^T] Q)
          let cross = 0, tr = 0; for (let a = 0; a < K; a++) cross += v.W[j * K + a] * r[a];
          for (let t = 0; t < KK; t++) tr += v.WW[j * KK + t] * Q[t];
          const e = Math.max(yy - 2 * cross + tr, 1e-10);
          v.tau[j] = (a0 + c / 2) / (b0 + e / 2); sse += e; sobs += c;
        }
        for (let a = 0; a < K; a++) { let s = 0; for (let j = 0; j < v.k; j++) s += v.WW[j * KK + a * K + a]; v.alpha[a] = (a0 + v.k / 2) / (b0 + s / 2); }
      }
      let dz = 0, zn = 0; for (let t = 0; t < n * K; t++) { dz += (Z[t] - Zprev[t]) ** 2; zn += Z[t] * Z[t]; }
      const crit = Math.sqrt(dz / Math.max(zn, 1e-12)); elboTrace.push(sse / sobs);
      if (opts.onIter && iter % 10 === 0) opts.onIter(iter, crit);
      if (crit < tol && iter > 10) { converged = true; break; }
    }
    // variance explained (MOFA2 definition, observed entries, centered data)
    const r2 = V.map(v => {
      let tss = 0; for (let i = 0; i < n * v.k; i++) if (v.obs[i]) tss += v.M[i] * v.M[i];
      const perFactor = []; let rssAll = 0;
      for (let a = 0; a < K; a++) { let rss = 0; for (let i = 0; i < n; i++) for (let j = 0; j < v.k; j++) { const idx = i * v.k + j; if (!v.obs[idx]) continue; rss += (v.M[idx] - Z[i * K + a] * v.W[j * K + a]) ** 2; } perFactor.push(Math.max(0, 1 - rss / tss)); }
      for (let i = 0; i < n; i++) for (let j = 0; j < v.k; j++) { const idx = i * v.k + j; if (!v.obs[idx]) continue; let f = 0; for (let a = 0; a < K; a++) f += Z[i * K + a] * v.W[j * K + a]; rssAll += (v.M[idx] - f) ** 2; }
      return { perFactor, total: Math.max(0, 1 - rssAll / tss) };
    });
    // order factors by summed R2, fix sign (largest-|weight| feature in the dominant view positive)
    const order = Array.from({ length: K }, (_, a) => a).sort((a, b) => r2.reduce((s, x) => s + x.perFactor[b], 0) - r2.reduce((s, x) => s + x.perFactor[a], 0));
    const signs = order.map(a => {
      let best = 0, bestV = 0; V.forEach((v, m) => { if (r2[m].perFactor[a] > bestV) { bestV = r2[m].perFactor[a]; best = m; } });
      const v = V[best]; let mj = 0; for (let j = 0; j < v.k; j++) if (Math.abs(v.W[j * K + a]) > Math.abs(v.W[mj * K + a])) mj = j;
      return v.W[mj * K + a] < 0 ? -1 : 1;
    });
    const Zo = Array.from({ length: n }, (_, i) => order.map((a, o) => Z[i * K + a] * signs[o]));
    const Wo = V.map(v => Array.from({ length: v.k }, (_, j) => order.map((a, o) => v.W[j * K + a] * signs[o])));
    return {
      Z: Zo, W: Wo, K, iterations: iter + 1, converged, trace: elboTrace,
      r2: r2.map(x => ({ perFactor: order.map(a => x.perFactor[a]), total: x.total })),
      alpha: V.map(v => order.map(a => v.alpha[a]))
    };
  }

  // match factors between two fits by |corr| (greedy)
  function matchFactors(ZA, ZB) {
    const KA = ZA[0].length, KB = ZB[0].length, n = ZA.length, C = [];
    for (let a = 0; a < KA; a++) { C.push([]); for (let b = 0; b < KB; b++) C[a].push(pearson(ZA.map(r => r[a]), ZB.map(r => r[b])).r); }
    const used = new Set(), out = new Array(KA).fill(null);
    const pairs = []; for (let a = 0; a < KA; a++) for (let b = 0; b < KB; b++) pairs.push([Math.abs(C[a][b]), a, b]);
    pairs.sort((x, y) => y[0] - x[0]); const usedA = new Set();
    for (const [v, a, b] of pairs) { if (usedA.has(a) || used.has(b)) continue; usedA.add(a); used.add(b); out[a] = { b, r: C[a][b] }; }
    return { matrix: C, match: out };
  }

  // ---------- association of a score vector with covariates ----------
  function associate(score, cov) { // cov: {numeric, codes}
    const y = [], g = []; for (let i = 0; i < score.length; i++) { const c = cov.codes[i]; if (c == null || c !== c && cov.numeric) continue; y.push(score[i]); g.push(c); }
    if (cov.numeric) { const s = spearman(y, g); return { stat: s.r, p: s.p, kind: 'Spearman ρ' }; }
    const lv = [...new Set(g)]; if (lv.length < 2) return { stat: NaN, p: NaN, kind: 'ANOVA' };
    const gm = mean(y); let ssb = 0, ssw = 0;
    lv.forEach(l => { const v = y.filter((_, i) => g[i] === l); const m = mean(v); ssb += v.length * (m - gm) ** 2; v.forEach(x => ssw += (x - m) ** 2); });
    const d1 = lv.length - 1, d2 = y.length - lv.length; const F = (ssb / d1) / (ssw / d2);
    return { stat: ssb / (ssb + ssw), p: pf(F, d1, d2), kind: 'η²' };
  }

  // ---------- enrichment ----------
  function ora(selected, universe, sets, minSize, maxSize) {
    const U = new Set(universe), S = new Set(selected.filter(g => U.has(g))); const N = U.size, n = S.size; const out = [];
    for (const [id, set] of Object.entries(sets)) {
      const members = set.genes.filter(g => U.has(g)); if (members.length < minSize || members.length > maxSize) continue;
      const hits = members.filter(g => S.has(g)); if (hits.length < 2) continue;
      out.push({ id, name: set.name, size: members.length, hits: hits.length, genes: hits, ratio: (hits.length / n) / (members.length / N), p: phyperUpper(hits.length, members.length, N, n) });
    }
    const q = bh(out.map(r => r.p)); out.forEach((r, i) => r.q = q[i]); out.sort((a, b) => a.p - b.p);
    return { results: out, nSelected: n, nUniverse: N };
  }
  // rank-based: two-sided Mann-Whitney of set members' scores vs the rest; effect = AUC-0.5 signed
  function rankTest(scores, sets, minSize, maxSize) { // scores: {gene: value}
    const genes = Object.keys(scores), vals = genes.map(g => scores[g]), r = ranks(vals), N = genes.length; const pos = new Map(genes.map((g, i) => [g, i]));
    // tie correction term
    const counts = new Map(); vals.forEach(v => counts.set(v, (counts.get(v) || 0) + 1)); let tie = 0; counts.forEach(c => { if (c > 1) tie += c * c * c - c; });
    const out = [];
    for (const [id, set] of Object.entries(sets)) {
      const idx = set.genes.filter(g => pos.has(g)).map(g => pos.get(g)); const n1 = idx.length; if (n1 < minSize || n1 > maxSize) continue;
      const n2 = N - n1; let R1 = 0; idx.forEach(i => R1 += r[i]);
      const U = R1 - n1 * (n1 + 1) / 2, mu = n1 * n2 / 2, sd = Math.sqrt(n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1))));
      const z = (U - mu) / sd; const auc = U / (n1 * n2);
      const lead = idx.map(i => [genes[i], vals[i]]).sort((a, b) => (z > 0 ? b[1] - a[1] : a[1] - b[1])).slice(0, 6).map(x => x[0]);
      out.push({ id, name: set.name, size: n1, z, auc, p: 2 * pnormUpper(Math.abs(z)), genes: lead });
    }
    const q = bh(out.map(x => x.p)); out.forEach((x, i) => x.q = q[i]); out.sort((a, b) => a.p - b.p);
    return { results: out, nRanked: N };
  }

  // ---------- classification with cross-validation ----------
  function logisticRidge(X, y, n, p, lambda) { // X n x p (standardized), y 0/1 -> beta (p+1)
    const q = p + 1; let beta = new Float64Array(q);
    for (let it = 0; it < 30; it++) {
      const H = new Float64Array(q * q), g = new Float64Array(q);
      for (let i = 0; i < n; i++) {
        let eta = beta[0]; for (let j = 0; j < p; j++) eta += beta[j + 1] * X[i * p + j];
        const mu = 1 / (1 + Math.exp(-eta)), w = Math.max(mu * (1 - mu), 1e-6), e = y[i] - mu;
        const row = [1]; for (let j = 0; j < p; j++) row.push(X[i * p + j]);
        for (let a = 0; a < q; a++) { g[a] += row[a] * e; for (let b = 0; b <= a; b++) H[a * q + b] += w * row[a] * row[b]; }
      }
      for (let a = 0; a < q; a++) { for (let b = a + 1; b < q; b++) H[a * q + b] = H[b * q + a]; if (a > 0) { H[a * q + a] += lambda; g[a] -= lambda * beta[a]; } else H[0] += 1e-8; }
      const step = solveSPD(H, g, q); let mx = 0; for (let a = 0; a < q; a++) { beta[a] += step[a]; mx = Math.max(mx, Math.abs(step[a])); }
      if (mx < 1e-7) break;
    }
    return beta;
  }
  function auc(scores, labels) {
    const pos = [], neg = []; scores.forEach((s, i) => (labels[i] ? pos : neg).push(s)); if (!pos.length || !neg.length) return NaN;
    const all = scores.slice(), r = ranks(all); let R = 0; labels.forEach((l, i) => { if (l) R += r[i]; });
    return (R - pos.length * (pos.length + 1) / 2) / (pos.length * neg.length);
  }
  function rocCurve(scores, labels) {
    const idx = scores.map((s, i) => i).sort((a, b) => scores[b] - scores[a]); const P = labels.filter(Boolean).length, N = labels.length - P;
    const pts = [[0, 0]]; let tp = 0, fp = 0;
    for (let k = 0; k < idx.length; k++) { if (labels[idx[k]]) tp++; else fp++; if (k === idx.length - 1 || scores[idx[k + 1]] !== scores[idx[k]]) pts.push([fp / N, tp / P]); }
    return pts;
  }
  /* blocks: [{name, M (n x k, NaN allowed)}], y: 0/1; options: {folds, repeats, topN, leaky, lambda, seed, permute} */
  function crossValidate(blocks, y, n, opt) {
    const rand = rng(opt.seed || 557); const out = blocks.map(() => ({ scores: new Float64Array(n).fill(0), counts: new Float64Array(n).fill(0), aucs: [] }));
    let yy = y.slice(); if (opt.permute) { const p = shuffleIdx(n, rand); yy = p.map(i => y[i]); }
    // leaky: select features on all samples before CV
    const leakySel = blocks.map(b => opt.leaky ? topFeaturesByT(b.M, n, b.k, yy, range(n), opt.topN) : null);
    for (let rep = 0; rep < opt.repeats; rep++) {
      // stratified folds
      const i1 = shuffleIdx(n, rand).filter(i => yy[i]), i0 = shuffleIdx(n, rand).filter(i => !yy[i]); const fold = new Array(n);
      i1.forEach((i, k) => fold[i] = k % opt.folds); i0.forEach((i, k) => fold[i] = k % opt.folds);
      const repScores = blocks.map(() => new Float64Array(n));
      for (let f = 0; f < opt.folds; f++) {
        const tr = range(n).filter(i => fold[i] !== f), te = range(n).filter(i => fold[i] === f);
        blocks.forEach((b, bi) => {
          const sel = opt.leaky ? leakySel[bi] : topFeaturesByT(b.M, n, b.k, yy, tr, opt.topN);
          const { Xtr, Xte } = standardizeSplit(b.M, n, b.k, sel, tr, te);
          const beta = logisticRidge(Xtr, tr.map(i => yy[i]), tr.length, sel.length, opt.lambda || 1);
          te.forEach((i, r) => { let eta = beta[0]; for (let j = 0; j < sel.length; j++) eta += beta[j + 1] * Xte[r * sel.length + j]; repScores[bi][i] = eta; });
        });
      }
      blocks.forEach((b, bi) => { out[bi].aucs.push(auc(Array.from(repScores[bi]), yy)); for (let i = 0; i < n; i++) { out[bi].scores[i] += repScores[bi][i]; out[bi].counts[i]++; } });
    }
    return blocks.map((b, bi) => {
      const s = Array.from(out[bi].scores).map((v, i) => v / out[bi].counts[i]); const a = out[bi].aucs;
      const m = mean(a), sd = Math.sqrt(a.reduce((x, v) => x + (v - m) ** 2, 0) / Math.max(1, a.length - 1));
      return { name: b.name, auc: m, aucSD: sd, aucs: a, roc: rocCurve(s, yy), meanScore: s };
    });
  }
  function range(n) { return Array.from({ length: n }, (_, i) => i); }
  function topFeaturesByT(M, n, k, y, rows, topN) {
    const t = []; for (let j = 0; j < k; j++) {
      const a = [], b = []; rows.forEach(i => { const v = M[i * k + j]; if (v === v) (y[i] ? b : a).push(v); });
      if (a.length < 2 || b.length < 2) { t.push([j, 0]); continue; }
      const ma = mean(a), mb = mean(b), va = a.reduce((s, v) => s + (v - ma) ** 2, 0) / (a.length - 1), vb = b.reduce((s, v) => s + (v - mb) ** 2, 0) / (b.length - 1);
      t.push([j, Math.abs(mb - ma) / Math.sqrt(va / a.length + vb / b.length + 1e-12)]);
    }
    return t.sort((x, z) => z[1] - x[1]).slice(0, Math.min(topN, k)).map(x => x[0]);
  }
  function standardizeSplit(M, n, k, sel, tr, te) {
    const p = sel.length, Xtr = new Float64Array(tr.length * p), Xte = new Float64Array(te.length * p);
    sel.forEach((j, c) => {
      const v = tr.map(i => M[i * k + j]).filter(x => x === x); const mu = v.length ? mean(v) : 0;
      const sd = v.length > 1 ? Math.sqrt(v.reduce((s, x) => s + (x - mu) ** 2, 0) / (v.length - 1)) || 1 : 1;
      tr.forEach((i, r) => { const x = M[i * k + j]; Xtr[r * p + c] = x === x ? (x - mu) / sd : 0; });
      te.forEach((i, r) => { const x = M[i * k + j]; Xte[r * p + c] = x === x ? (x - mu) / sd : 0; });
    });
    return { Xtr, Xte };
  }

  const Engine = {
    rng, shuffleIdx, topEig, lgamma, digamma, trigamma, trigammaInverse, ibeta, pt2, pf, pnormUpper, phyperUpper, bh,
    mean, median, quantile, ranks, pearson, spearman, corP, eigSym, cholInv,
    toMatrix, column, subsetCols, subsetRows, colStats, prepareView, imputeColumns, pca, fitFDist, featureTests,
    variancePartition, gfa, matchFactors, associate, ora, rankTest, crossValidate, auc, rocCurve, range
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = Engine; else root.Engine = Engine;
})(typeof self !== 'undefined' ? self : this);
