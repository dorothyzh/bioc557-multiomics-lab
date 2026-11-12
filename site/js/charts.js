/* Publication-style SVG charts (d3 v7).
 * Every chart is drawn in print units: 1 user unit = 1 pt. Width presets follow journal columns
 * (single 89 mm, 1.5-column 120 mm, double 183 mm). On screen the SVG scales to its container;
 * exported SVG/PNG keep the physical size, so fonts stay at 6–7 pt in print.
 */
(function (root) {
  'use strict';
  const d3 = root.d3;
  const MM = 72 / 25.4;
  const WIDTHS = { single: 89 * MM, onehalf: 120 * MM, double: 183 * MM };
  const FONT = '"Helvetica Neue", Helvetica, Arial, sans-serif';
  const INK = '#1a1a1a', INK2 = '#4d4d4d', MUTED = '#8c8c8c', RULE = '#bdbdbd', GRID = '#ececec', NS = '#c4c4c4';
  const FS = { tick: 6.2, label: 7, legend: 6.4, note: 6, annot: 6 };
  const LW = { axis: 0.6, tick: 0.6, line: 1.0, thin: 0.5 };
  const tip = (() => {
    let el = null;
    return {
      show(html, ev) {
        if (!el) { el = document.createElement('div'); el.className = 'viz-tip'; document.body.appendChild(el); }
        el.innerHTML = html; el.style.opacity = 1;
        const pad = 14, r = el.getBoundingClientRect();
        let x = ev.clientX + pad, y = ev.clientY + pad;
        if (x + r.width > window.innerWidth - 8) x = ev.clientX - r.width - pad;
        if (y + r.height > window.innerHeight - 8) y = ev.clientY - r.height - pad;
        el.style.transform = `translate(${x}px, ${y}px)`;
      },
      hide() { if (el) el.style.opacity = 0; }
    };
  })();
  const fmt = (v, d = 2) => (v == null || !isFinite(v) ? '—' : Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(1) : (+v).toFixed(d));
  const fmtP = v => (v == null || !isFinite(v) ? '—' : v < 1e-4 ? v.toExponential(1) : v.toPrecision(2));
  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

  function frame(el, o) {
    const W = typeof o.width === 'number' ? o.width : WIDTHS[o.width || 'single'];
    const H = o.height || W * (o.aspect || 0.78);
    const m = Object.assign({ t: 8, r: 8, b: 30, l: 36 }, o.margin || {});
    d3.select(el).selectAll('*').remove();
    const svg = d3.select(el).append('svg')
      .attr('xmlns', 'http://www.w3.org/2000/svg').attr('viewBox', `0 0 ${W} ${H}`)
      .attr('data-width-mm', (W / MM).toFixed(1)).attr('data-height-mm', (H / MM).toFixed(1))
      .attr('font-family', FONT).attr('font-size', FS.tick).attr('fill', INK)
      .style('width', '100%').style('max-width', `${Math.round(W * 2.05)}px`).style('height', 'auto').style('display', 'block').style('margin', '0 auto');
    svg.append('rect').attr('class', 'bg').attr('width', W).attr('height', H).attr('fill', '#ffffff');
    const g = svg.append('g').attr('transform', `translate(${m.l},${m.t})`);
    return { svg, g, W, H, m, iw: W - m.l - m.r, ih: H - m.t - m.b };
  }
  function styleAxis(sel) {
    sel.attr('font-family', FONT).attr('font-size', FS.tick);
    sel.selectAll('.domain').attr('stroke', INK).attr('stroke-width', LW.axis);
    sel.selectAll('.tick line').attr('stroke', INK).attr('stroke-width', LW.tick);
    sel.selectAll('.tick text').attr('fill', INK);
    return sel;
  }
  function axes(f, x, y, o) {
    o = o || {};
    const nx = o.xTicks || Math.max(3, Math.round(f.iw / 55)), ny = o.yTicks || Math.max(3, Math.round(f.ih / 40));
    if (o.grid) {
      f.g.append('g').attr('class', 'grid').selectAll('line').data(y.ticks ? y.ticks(ny) : []).join('line')
        .attr('x1', 0).attr('x2', f.iw).attr('y1', d => y(d)).attr('y2', d => y(d)).attr('stroke', GRID).attr('stroke-width', LW.thin);
    }
    const xa = f.g.append('g').attr('transform', `translate(0,${f.ih + (o.offset || 0)})`).call(o.xAxis || d3.axisBottom(x).ticks(nx, o.xFormat).tickSizeOuter(0).tickSize(2.5).tickPadding(2));
    const ya = f.g.append('g').attr('transform', `translate(${-(o.offset || 0)},0)`).call(o.yAxis || d3.axisLeft(y).ticks(ny, o.yFormat).tickSizeOuter(0).tickSize(2.5).tickPadding(2));
    styleAxis(xa); styleAxis(ya);
    if (o.xLabel) f.g.append('text').attr('x', f.iw / 2).attr('y', f.ih + (o.offset || 0) + 19).attr('text-anchor', 'middle').attr('font-size', FS.label).text(o.xLabel);
    if (o.yLabel) f.g.append('text').attr('transform', `translate(${-(o.offset || 0) - (o.yLabelOffset || 25)},${f.ih / 2}) rotate(-90)`).attr('text-anchor', 'middle').attr('font-size', FS.label).text(o.yLabel);
    return { xa, ya };
  }
  // legend in a row (top-right inside plotting area by default)
  function legend(f, items, o) {
    o = o || {}; const g = f.g.append('g').attr('class', 'legend');
    let x = 0; const parts = [];
    items.forEach(it => {
      const gi = g.append('g').attr('transform', `translate(${x},0)`);
      if (it.type === 'line') gi.append('line').attr('x1', 0).attr('x2', 9).attr('y1', 0).attr('y2', 0).attr('stroke', it.color).attr('stroke-width', 1.2).attr('stroke-dasharray', it.dash || null);
      else if (it.type === 'square') gi.append('rect').attr('x', 0).attr('y', -3).attr('width', 6).attr('height', 6).attr('rx', 0.8).attr('fill', it.color);
      else symbol(gi, it.shape || 'circle', 3.5, 0, 2.3).attr('fill', it.color).attr('stroke', '#fff').attr('stroke-width', 0.4);
      const tx = gi.append('text').attr('x', it.type === 'line' ? 12 : 9).attr('y', 2.2).attr('font-size', FS.legend).attr('fill', INK).text(it.label);
      const w = (it.type === 'line' ? 12 : 9) + textWidth(it.label, FS.legend) + 8; x += w; parts.push(w);
    });
    const total = x - 8;
    let ax = o.align === 'left' ? (o.x || 0) : (o.x != null ? o.x : f.iw - total);
    if (ax < -f.m.l + 2) ax = -f.m.l + 2;
    g.attr('transform', `translate(${ax},${o.y != null ? o.y : -2})`);
    return g;
  }
  function textWidth(s, size) { return String(s).length * size * 0.54; }
  function symbol(sel, shape, cx, cy, r) {
    if (shape === 'triangle') return sel.append('path').attr('d', `M${cx},${cy - r * 1.15} L${cx + r * 1.1},${cy + r * 0.75} L${cx - r * 1.1},${cy + r * 0.75} Z`);
    if (shape === 'square') return sel.append('rect').attr('x', cx - r * 0.9).attr('y', cy - r * 0.9).attr('width', r * 1.8).attr('height', r * 1.8);
    if (shape === 'diamond') return sel.append('path').attr('d', `M${cx},${cy - r * 1.2} L${cx + r * 1.1},${cy} L${cx},${cy + r * 1.2} L${cx - r * 1.1},${cy} Z`);
    return sel.append('circle').attr('cx', cx).attr('cy', cy).attr('r', r);
  }
  function symbolPath(shape, r) {
    const map = { circle: d3.symbolCircle, triangle: d3.symbolTriangle, square: d3.symbolSquare, diamond: d3.symbolDiamond };
    const area = { circle: Math.PI * r * r, triangle: r * r * 2.1, square: r * r * 3.2, diamond: r * r * 2.6 }[shape || 'circle'];
    return d3.symbol(map[shape || 'circle'] || d3.symbolCircle, area)();
  }
  function pad(ext, f) { const span = (ext[1] - ext[0]) || 1; return [ext[0] - span * (f || 0.06), ext[1] + span * (f || 0.06)]; }

  // 95% normal-theory ellipse from points
  function ellipsePath(pts, x, y) {
    if (pts.length < 4) return null;
    const mx = d3.mean(pts, p => p.x), my = d3.mean(pts, p => p.y);
    let sxx = 0, syy = 0, sxy = 0; pts.forEach(p => { sxx += (p.x - mx) ** 2; syy += (p.y - my) ** 2; sxy += (p.x - mx) * (p.y - my); });
    const n1 = pts.length - 1; sxx /= n1; syy /= n1; sxy /= n1;
    const tr = sxx + syy, det = sxx * syy - sxy * sxy, l1 = tr / 2 + Math.sqrt(Math.max(0, tr * tr / 4 - det)), l2 = tr / 2 - Math.sqrt(Math.max(0, tr * tr / 4 - det));
    const ang = Math.atan2(l1 - sxx, sxy || 1e-12), c = Math.sqrt(5.991);
    const out = []; for (let t = 0; t <= 64; t++) { const th = (t / 64) * 2 * Math.PI; const ex = c * Math.sqrt(Math.max(l1, 0)) * Math.cos(th), ey = c * Math.sqrt(Math.max(l2, 0)) * Math.sin(th); out.push([x(mx + ex * Math.cos(ang) - ey * Math.sin(ang)), y(my + ex * Math.sin(ang) + ey * Math.cos(ang))]); }
    return d3.line()(out) + 'Z';
  }

  // greedy label placement avoiding overlaps
  function placeLabels(g, items, bounds, size) {
    size = size || FS.annot; const boxes = [];
    const overlaps = b => boxes.some(o => !(b.x2 < o.x1 || b.x1 > o.x2 || b.y2 < o.y1 || b.y1 > o.y2));
    items.forEach(it => {
      const w = textWidth(it.text, size), h = size * 1.05;
      const cands = [[5, -3, 'start'], [-5, -3, 'end'], [5, 7, 'start'], [-5, 7, 'end'], [0, -6, 'middle'], [0, 10, 'middle'], [9, -9, 'start'], [-9, -9, 'end'], [9, 12, 'start'], [-9, 12, 'end']];
      for (const [dx, dy, anchor] of cands) {
        const lx = it.x + dx, ly = it.y + dy;
        const x1 = anchor === 'start' ? lx : anchor === 'end' ? lx - w : lx - w / 2;
        const b = { x1, x2: x1 + w, y1: ly - h * 0.8, y2: ly + h * 0.2 };
        if (b.x1 < bounds[0] || b.x2 > bounds[2] || b.y1 < bounds[1] || b.y2 > bounds[3] || overlaps(b)) continue;
        boxes.push(b);
        if (Math.abs(dx) > 6 || Math.abs(dy) > 8) g.append('line').attr('x1', it.x).attr('y1', it.y).attr('x2', anchor === 'start' ? x1 - 0.5 : anchor === 'end' ? x1 + w + 0.5 : lx).attr('y2', ly - h * 0.3).attr('stroke', MUTED).attr('stroke-width', 0.35);
        g.append('text').attr('x', lx).attr('y', ly).attr('text-anchor', anchor).attr('font-size', size).attr('fill', INK).attr('font-style', it.italic ? 'italic' : null).text(it.text);
        break;
      }
    });
  }

  // ---------- scatter ----------
  function scatter(el, o) {
    const f = frame(el, { width: o.width, aspect: o.aspect || 0.8, margin: Object.assign({ t: 14, l: 34, b: 28 }, o.margin) });
    const pts = o.points.filter(p => isFinite(p.x) && isFinite(p.y));
    let xd = pad(d3.extent(pts, p => p.x)), yd = pad(d3.extent(pts, p => p.y));
    if (o.square) { const lo = Math.min(xd[0], yd[0]), hi = Math.max(xd[1], yd[1]); xd = [lo, hi]; yd = [lo, hi]; }
    const x = d3.scaleLinear().domain(xd).range([0, f.iw]).nice(), y = d3.scaleLinear().domain(yd).range([f.ih, 0]).nice();
    axes(f, x, y, { xLabel: o.xLabel, yLabel: o.yLabel, yLabelOffset: o.yLabelOffset });
    if (o.zeroLines) {
      if (x.domain()[0] < 0 && x.domain()[1] > 0) f.g.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 0).attr('y2', f.ih).attr('stroke', RULE).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
      if (y.domain()[0] < 0 && y.domain()[1] > 0) f.g.append('line').attr('x1', 0).attr('x2', f.iw).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', RULE).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    }
    if (o.diagonal) f.g.append('line').attr('x1', x(Math.max(x.domain()[0], y.domain()[0]))).attr('y1', y(Math.max(x.domain()[0], y.domain()[0]))).attr('x2', x(Math.min(x.domain()[1], y.domain()[1]))).attr('y2', y(Math.min(x.domain()[1], y.domain()[1]))).attr('stroke', RULE).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    const clip = f.g.append('g');
    if (o.ellipses) o.groups.forEach(gr => { const pth = ellipsePath(pts.filter(p => p.group === gr.name), x, y); if (pth) clip.append('path').attr('d', pth).attr('fill', gr.color).attr('fill-opacity', 0.07).attr('stroke', gr.color).attr('stroke-width', 0.6).attr('stroke-opacity', 0.8); });
    if (o.fit) { // least squares line per all points
      const n = pts.length, mx = d3.mean(pts, p => p.x), my = d3.mean(pts, p => p.y); let sxy = 0, sxx = 0; pts.forEach(p => { sxy += (p.x - mx) * (p.y - my); sxx += (p.x - mx) ** 2; });
      const b = sxy / sxx, a = my - b * mx, [x0, x1] = d3.extent(pts, p => p.x);
      clip.append('line').attr('x1', x(x0)).attr('y1', y(a + b * x0)).attr('x2', x(x1)).attr('y2', y(a + b * x1)).attr('stroke', INK2).attr('stroke-width', 0.8);
    }
    const color = gname => (o.groups.find(g => g.name === gname) || { color: MUTED }).color;
    const shape = gname => (o.groups.find(g => g.name === gname) || {}).shape || 'circle';
    const r = o.r || (pts.length > 150 ? 1.7 : 2.3);
    clip.selectAll('path.pt').data(pts).join('path').attr('class', 'pt')
      .attr('d', p => symbolPath(p.shape || shape(p.group), r)).attr('transform', p => `translate(${x(p.x)},${y(p.y)})`)
      .attr('fill', p => p.color || color(p.group)).attr('fill-opacity', o.opacity || 0.9).attr('stroke', '#fff').attr('stroke-width', 0.4)
      .on('mousemove', (ev, p) => tip.show(o.tooltip ? o.tooltip(p) : `<b>${esc(p.label || p.id)}</b><br>${esc(p.group || '')}<br>${fmt(p.x)}, ${fmt(p.y)}`, ev))
      .on('mouseleave', tip.hide).on('click', (ev, p) => o.onClick && o.onClick(p)).style('cursor', o.onClick ? 'pointer' : null);
    if (o.labels && o.labels.length) placeLabels(f.g, o.labels.map(l => ({ x: x(l.x), y: y(l.y), text: l.text, italic: l.italic })), [0, 0, f.iw, f.ih]);
    if (o.groups && o.groups.length > 1 && o.legend !== false) legend(f, o.groups.map(g => ({ label: g.label || g.name, color: g.color, shape: g.shape })), { y: -7 });
    if (o.stat) f.g.append('text').attr('x', 2).attr('y', 4).attr('font-size', FS.note).attr('fill', INK2).text(o.stat);
    return f;
  }

  // ---------- volcano ----------
  function volcano(el, o) {
    const f = frame(el, { width: o.width, aspect: o.aspect || 0.8, margin: { t: 14, l: 30, b: 28, r: 8 } });
    const pts = o.points.filter(p => isFinite(p.x) && isFinite(p.y));
    const xm = d3.max(pts, p => Math.abs(p.x)) || 1;
    const x = d3.scaleLinear().domain([-xm * 1.08, xm * 1.08]).range([0, f.iw]).nice(), y = d3.scaleLinear().domain([0, (d3.max(pts, p => p.y) || 1) * 1.06]).range([f.ih, 0]).nice();
    axes(f, x, y, { xLabel: o.xLabel, yLabel: o.yLabel || '−log₁₀ P', yLabelOffset: 20 });
    if (o.pLine != null) f.g.append('line').attr('x1', 0).attr('x2', f.iw).attr('y1', y(o.pLine)).attr('y2', y(o.pLine)).attr('stroke', MUTED).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    if (o.fcLine) [-o.fcLine, o.fcLine].forEach(v => f.g.append('line').attr('x1', x(v)).attr('x2', x(v)).attr('y1', 0).attr('y2', f.ih).attr('stroke', MUTED).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5'));
    const col = p => (p.sig ? (p.x > 0 ? o.upColor : o.downColor) : NS);
    const sorted = pts.slice().sort((a, b) => (a.sig ? 1 : 0) - (b.sig ? 1 : 0));
    f.g.selectAll('circle').data(sorted).join('circle').attr('cx', p => x(p.x)).attr('cy', p => y(p.y)).attr('r', p => (p.sig ? 1.9 : 1.4))
      .attr('fill', col).attr('fill-opacity', p => (p.sig ? 0.9 : 0.7)).attr('stroke', '#fff').attr('stroke-width', 0.3)
      .on('mousemove', (ev, p) => tip.show(o.tooltip(p), ev)).on('mouseleave', tip.hide).on('click', (ev, p) => o.onClick && o.onClick(p)).style('cursor', o.onClick ? 'pointer' : null);
    if (o.labels) placeLabels(f.g, o.labels.map(p => ({ x: x(p.x), y: y(p.y), text: p.label, italic: o.italic })), [0, 0, f.iw, f.ih]);
    const nu = pts.filter(p => p.sig && p.x > 0).length, nd = pts.filter(p => p.sig && p.x < 0).length;
    f.g.append('text').attr('x', f.iw - 2).attr('y', 3).attr('text-anchor', 'end').attr('font-size', FS.note).attr('fill', o.upColor).text(`${nu} higher in ${o.upLabel}`);
    f.g.append('text').attr('x', 2).attr('y', 3).attr('font-size', FS.note).attr('fill', o.downColor).text(`${nd} higher in ${o.downLabel}`);
    return f;
  }

  // ---------- bars per sample (vertical, thin) ----------
  function sampleBars(el, o) {
    const f = frame(el, { width: o.width, aspect: o.aspect || 0.62, margin: { t: 14, l: 30, b: o.showNames ? 44 : 22, r: 6 } });
    const x = d3.scaleBand().domain(o.items.map(d => d.id)).range([0, f.iw]).paddingInner(o.items.length > 60 ? 0.12 : 0.2);
    const y = d3.scaleLinear().domain([0, o.max != null ? o.max : Math.max(o.minMax || 0, (d3.max(o.items, d => d.value) || 0) * 1.1)]).range([f.ih, 0]).nice();
    const xAxis = d3.axisBottom(x).tickSize(0).tickFormat(o.showNames ? null : '');
    axes(f, x, y, { xAxis, yFormat: o.yFormat, xLabel: o.showNames ? null : o.xLabel, yLabel: o.yLabel, grid: true, yLabelOffset: 22 });
    if (o.showNames) f.g.selectAll('.tick text').filter(function () { return this.parentNode.parentNode.getAttribute('transform')?.includes(`,${f.ih}`); })
      .attr('transform', 'rotate(-60)').attr('text-anchor', 'end').attr('dx', -2).attr('dy', 2).attr('font-size', 4.6);
    f.g.selectAll('rect.bar').data(o.items).join('rect').attr('class', 'bar').attr('x', d => x(d.id)).attr('width', x.bandwidth())
      .attr('y', d => y(Math.max(0, d.value))).attr('height', d => Math.abs(y(0) - y(d.value))).attr('fill', d => d.color).attr('rx', Math.min(0.8, x.bandwidth() / 3))
      .on('mousemove', (ev, d) => tip.show(o.tooltip(d), ev)).on('mouseleave', tip.hide);
    if (o.legendItems) legend(f, o.legendItems.map(l => ({ ...l, type: 'square' })), { y: -7 });
    if (o.hline != null) f.g.append('line').attr('x1', 0).attr('x2', f.iw).attr('y1', y(o.hline)).attr('y2', y(o.hline)).attr('stroke', INK2).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    return f;
  }

  // ---------- per-sample box summaries ----------
  function sampleBoxes(el, o) {
    const f = frame(el, { width: o.width, aspect: o.aspect || 0.62, margin: { t: 14, l: 30, b: 22, r: 6 } });
    const x = d3.scaleBand().domain(o.items.map(d => d.id)).range([0, f.iw]).paddingInner(0.25);
    const lo = d3.min(o.items, d => d.q[0]), hi = d3.max(o.items, d => d.q[4]);
    const y = d3.scaleLinear().domain(pad([lo, hi], 0.04)).range([f.ih, 0]).nice();
    axes(f, x, y, { xAxis: d3.axisBottom(x).tickSize(0).tickFormat(''), xLabel: o.xLabel, yLabel: o.yLabel, grid: true, yLabelOffset: 22 });
    const bw = x.bandwidth(), gs = f.g.selectAll('g.box').data(o.items).join('g').attr('class', 'box').attr('transform', d => `translate(${x(d.id)},0)`)
      .on('mousemove', (ev, d) => tip.show(o.tooltip(d), ev)).on('mouseleave', tip.hide);
    gs.append('line').attr('x1', bw / 2).attr('x2', bw / 2).attr('y1', d => y(d.q[0])).attr('y2', d => y(d.q[4])).attr('stroke', d => d.color).attr('stroke-width', Math.max(0.3, Math.min(0.6, bw / 4)));
    gs.append('rect').attr('x', 0).attr('width', bw).attr('y', d => y(d.q[3])).attr('height', d => Math.max(0.5, y(d.q[1]) - y(d.q[3]))).attr('fill', d => d.color).attr('fill-opacity', 0.85).attr('rx', Math.min(0.6, bw / 4));
    gs.append('line').attr('x1', 0).attr('x2', bw).attr('y1', d => y(d.q[2])).attr('y2', d => y(d.q[2])).attr('stroke', '#fff').attr('stroke-width', 0.6);
    if (o.legendItems) legend(f, o.legendItems.map(l => ({ ...l, type: 'square' })), { y: -7 });
    return f;
  }

  // ---------- histogram ----------
  function histogram(el, o) {
    const f = frame(el, { width: o.width, aspect: o.aspect || 0.66, margin: { t: 14, l: 30, b: 28, r: 8 } });
    const vals = o.values.filter(isFinite), dom = o.domain || d3.extent(vals);
    const x = d3.scaleLinear().domain(dom).range([0, f.iw]).nice();
    const bins = d3.bin().domain(x.domain()).thresholds(x.ticks(o.bins || 30))(vals);
    const y = d3.scaleLinear().domain([0, d3.max(bins, b => b.length) * 1.12]).range([f.ih, 0]).nice();
    axes(f, x, y, { xLabel: o.xLabel, yLabel: o.yLabel || 'Count', yLabelOffset: 22 });
    f.g.selectAll('rect').data(bins).join('rect').attr('x', b => x(b.x0) + 0.25).attr('width', b => Math.max(0, x(b.x1) - x(b.x0) - 0.5))
      .attr('y', b => y(b.length)).attr('height', b => y(0) - y(b.length)).attr('fill', o.color || '#5b6b7a').attr('rx', 0.5)
      .on('mousemove', (ev, b) => tip.show(`${fmt(b.x0)} to ${fmt(b.x1)}<br><b>${b.length}</b> ${o.unit || 'items'}`, ev)).on('mouseleave', tip.hide);
    (o.markers || []).forEach(mk => {
      f.g.append('line').attr('x1', x(mk.value)).attr('x2', x(mk.value)).attr('y1', 0).attr('y2', f.ih).attr('stroke', mk.color || INK).attr('stroke-width', 0.8).attr('stroke-dasharray', mk.dash || null);
      f.g.append('text').attr('x', x(mk.value) + 2).attr('y', 4 + (mk.dy || 0)).attr('font-size', FS.note).attr('fill', INK).text(mk.label);
    });
    return f;
  }

  // ---------- horizontal lollipop (weights / loadings) ----------
  function lollipop(el, o) {
    const n = o.items.length, rowH = 9.2;
    const f = frame(el, { width: o.width, height: 26 + n * rowH + 22, margin: { t: 10, l: o.labelWidth || 70, b: 24, r: 10 } });
    const m = d3.max(o.items, d => Math.abs(d.value)) || 1;
    const x = d3.scaleLinear().domain(o.domain || [Math.min(0, -m * (o.symmetric === false ? 0 : 1)) * 1.08, m * 1.08]).range([0, f.iw]).nice();
    const y = d3.scaleBand().domain(o.items.map(d => d.label)).range([0, f.ih]).paddingInner(0.3);
    const xa = f.g.append('g').attr('transform', `translate(0,${f.ih})`).call(d3.axisBottom(x).ticks(4).tickSize(2.5).tickPadding(2).tickSizeOuter(0)); styleAxis(xa);
    f.g.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 0).attr('y2', f.ih).attr('stroke', INK).attr('stroke-width', LW.axis);
    f.g.append('text').attr('x', f.iw / 2).attr('y', f.ih + 19).attr('text-anchor', 'middle').attr('font-size', FS.label).text(o.xLabel || 'Weight');
    const rows = f.g.selectAll('g.row').data(o.items).join('g').attr('class', 'row').attr('transform', d => `translate(0,${y(d.label) + y.bandwidth() / 2})`)
      .on('mousemove', (ev, d) => tip.show(o.tooltip ? o.tooltip(d) : `<b>${esc(d.label)}</b><br>${fmt(d.value, 3)}`, ev)).on('mouseleave', tip.hide)
      .on('click', (ev, d) => o.onClick && o.onClick(d)).style('cursor', o.onClick ? 'pointer' : null);
    rows.append('line').attr('x1', x(0)).attr('x2', d => x(d.value)).attr('stroke', d => d.color || o.color).attr('stroke-width', 1.1);
    rows.append('circle').attr('cx', d => x(d.value)).attr('r', 2.2).attr('fill', d => d.color || o.color);
    rows.append('text').attr('x', -4).attr('y', 2.1).attr('text-anchor', 'end').attr('font-size', FS.tick).attr('font-style', o.italic ? 'italic' : null)
      .text(d => (d.label.length > (o.maxChars || 22) ? d.label.slice(0, (o.maxChars || 22) - 1) + '…' : d.label));
    return f;
  }

  // ---------- heatmap ----------
  function diverging(domainMax) {
    return d3.scaleLinear().domain([-domainMax, -domainMax / 2, 0, domainMax / 2, domainMax]).range(['#2166ac', '#92c5de', '#f7f7f7', '#f4a582', '#b2182b']).clamp(true);
  }
  function sequential(domain, hue) {
    const ramps = { blue: ['#f4f8fc', '#c6dbef', '#6baed6', '#2171b5', '#08306b'], teal: ['#f3faf7', '#c2e5d8', '#66bfa3', '#1b8a6b', '#0b4d3b'] };
    const r = ramps[hue || 'blue'];
    return d3.scaleLinear().domain(r.map((_, i) => domain[0] + (domain[1] - domain[0]) * i / (r.length - 1))).range(r).clamp(true);
  }
  function heatmap(el, o) {
    const nr = o.rows.length, nc = o.cols.length;
    const cell = o.cell || Math.min(18, Math.max(2.2, (WIDTHS[o.width || 'single'] - (o.ml || 60) - 46) / nc));
    const W = typeof o.width === 'number' ? o.width : WIDTHS[o.width || 'single'];
    const annH = (o.colAnnot ? 5 : 0), legW = 40;
    const ml = o.ml || 60, mt = (o.mt || 10) + annH, cw = Math.min(cell, (W - ml - legW - 8) / nc), ch = o.cellH || cw;
    const H = mt + nr * ch + (o.mb || 30);
    const f = frame(el, { width: W, height: H, margin: { t: mt, l: ml, b: o.mb || 30, r: W - ml - nc * cw } });
    const color = o.scale === 'diverging' ? diverging(o.max || 1) : sequential(o.domain || [0, o.max || 1], o.hue);
    const g = f.g;
    const cells = [];
    for (let i = 0; i < nr; i++) for (let j = 0; j < nc; j++) cells.push({ i, j, v: o.values[i][j] });
    g.selectAll('rect.cell').data(cells).join('rect').attr('class', 'cell').attr('x', d => d.j * cw).attr('y', d => d.i * ch).attr('width', cw + 0.02).attr('height', ch + 0.02)
      .attr('fill', d => (d.v == null || !isFinite(d.v) ? '#ffffff' : color(d.v)))
      .on('mousemove', (ev, d) => tip.show(o.tooltip ? o.tooltip(d) : `${esc(o.rows[d.i])} · ${esc(o.cols[d.j])}<br><b>${fmt(d.v)}</b>`, ev)).on('mouseleave', tip.hide);
    if (o.gaps) { g.selectAll('rect.cell').attr('stroke', '#fff').attr('stroke-width', o.gaps); }
    if (o.cellLabels) g.selectAll('text.cl').data(cells).join('text').attr('class', 'cl').attr('x', d => d.j * cw + cw / 2).attr('y', d => d.i * ch + ch / 2 + 2)
      .attr('text-anchor', 'middle').attr('font-size', Math.min(6.2, ch * 0.42)).attr('fill', d => (Math.abs(o.scale === 'diverging' ? d.v / (o.max || 1) : (d.v - (o.domain ? o.domain[0] : 0)) / ((o.domain ? o.domain[1] - o.domain[0] : o.max) || 1)) > 0.6 ? '#fff' : INK)).text(d => o.cellLabels(d.v, d));
    if (o.showRowLabels !== false && ch >= 4.5) g.selectAll('text.rl').data(o.rows).join('text').attr('class', 'rl').attr('x', -3).attr('y', (d, i) => i * ch + ch / 2 + 2).attr('text-anchor', 'end').attr('font-size', Math.min(FS.tick, ch * 0.85)).text(d => d);
    if (o.showColLabels !== false && cw >= 4.5) g.selectAll('text.cl2').data(o.cols).join('text').attr('class', 'cl2').attr('transform', (d, j) => `translate(${j * cw + cw / 2 + 2},${nr * ch + 3}) rotate(${o.colAngle != null ? o.colAngle : -45})`).attr('text-anchor', 'end').attr('font-size', Math.min(FS.tick, cw * 0.85)).text(d => d);
    if (o.colAnnot) { // array of {color}
      g.selectAll('rect.ca').data(o.colAnnot).join('rect').attr('x', (d, j) => j * cw).attr('y', -annH - 1).attr('width', cw + 0.02).attr('height', annH).attr('fill', d => d.color);
    }
    if (o.rowAnnot) g.selectAll('rect.ra').data(o.rowAnnot).join('rect').attr('x', -annH - 1).attr('y', (d, i) => i * ch).attr('width', annH).attr('height', ch + 0.02).attr('fill', d => d.color);
    // colour bar
    const lx = nc * cw + 10, lh = Math.min(60, nr * ch), lg = g.append('g').attr('transform', `translate(${lx},0)`);
    const id = 'grad' + Math.random().toString(36).slice(2, 8), defs = f.svg.append('defs'), grad = defs.append('linearGradient').attr('id', id).attr('x1', 0).attr('x2', 0).attr('y1', 1).attr('y2', 0);
    const dom = o.scale === 'diverging' ? [-(o.max || 1), o.max || 1] : (o.domain || [0, o.max || 1]);
    d3.range(0, 1.0001, 0.1).forEach(t => grad.append('stop').attr('offset', t).attr('stop-color', color(dom[0] + t * (dom[1] - dom[0]))));
    lg.append('rect').attr('width', 5).attr('height', lh).attr('fill', `url(#${id})`);
    const ls = d3.scaleLinear().domain(dom).range([lh, 0]);
    styleAxis(lg.append('g').attr('transform', 'translate(5,0)').call(d3.axisRight(ls).ticks(3).tickSize(2).tickPadding(1.5))).select('.domain').remove();
    if (o.legendTitle) lg.append('text').attr('x', 0).attr('y', lh + 9).attr('font-size', FS.note).text(o.legendTitle);
    return f;
  }

  // ---------- strip plot (values by group, deterministic jitter) ----------
  function strip(el, o) {
    const f = frame(el, { width: o.width, aspect: o.aspect || 0.75, margin: { t: 14, l: 30, b: 26, r: 6 } });
    const cats = o.categories; // [{name,label,color}]
    const x = d3.scaleBand().domain(cats.map(c => c.name)).range([0, f.iw]).paddingInner(0.35).paddingOuter(0.2);
    const vals = o.points.map(p => p.y).filter(isFinite);
    const y = d3.scaleLinear().domain(pad(d3.extent(vals), 0.08)).range([f.ih, 0]).nice();
    axes(f, x, y, { xAxis: d3.axisBottom(x).tickSize(0).tickPadding(4).tickFormat(c => (cats.find(k => k.name === c) || {}).label || c), yLabel: o.yLabel, yLabelOffset: 22 });
    const bw = x.bandwidth();
    cats.forEach(c => {
      const v = o.points.filter(p => p.cat === c.name).map(p => p.y).filter(isFinite).sort((a, b) => a - b); if (!v.length) return;
      const q1 = d3.quantile(v, 0.25), q2 = d3.quantile(v, 0.5), q3 = d3.quantile(v, 0.75), cx = x(c.name);
      f.g.append('rect').attr('x', cx + bw * 0.18).attr('width', bw * 0.64).attr('y', y(q3)).attr('height', Math.max(0.5, y(q1) - y(q3))).attr('fill', 'none').attr('stroke', INK2).attr('stroke-width', 0.6).attr('rx', 0.8);
      f.g.append('line').attr('x1', cx + bw * 0.18).attr('x2', cx + bw * 0.82).attr('y1', y(q2)).attr('y2', y(q2)).attr('stroke', INK).attr('stroke-width', 1);
    });
    const rand = (s => () => (s = (s * 16807) % 2147483647) / 2147483647)(7);
    const pts = o.points.map(p => ({ ...p, jx: (rand() - 0.5) * bw * 0.62 }));
    f.g.selectAll('circle').data(pts).join('circle').attr('cx', p => x(p.cat) + bw / 2 + p.jx).attr('cy', p => y(p.y)).attr('r', o.r || 1.9)
      .attr('fill', p => (cats.find(c => c.name === p.cat) || {}).color || MUTED).attr('fill-opacity', 0.75).attr('stroke', '#fff').attr('stroke-width', 0.35)
      .on('mousemove', (ev, p) => tip.show(o.tooltip ? o.tooltip(p) : `<b>${esc(p.label)}</b><br>${fmt(p.y)}`, ev)).on('mouseleave', tip.hide);
    if (o.stat) f.g.append('text').attr('x', f.iw).attr('y', 2).attr('text-anchor', 'end').attr('font-size', FS.note).attr('fill', INK2).text(o.stat);
    return f;
  }

  // ---------- grouped strip for many factors: heat-strip of scores ----------
  // ---------- line chart ----------
  function lines(el, o) {
    const f = frame(el, { width: o.width, aspect: o.aspect || 0.72, margin: { t: 16, l: 30, b: 28, r: 10 } });
    const allx = o.series.flatMap(s => s.points.map(p => p[0])), ally = o.series.flatMap(s => s.points.map(p => p[1])).concat(o.yInclude || []);
    const x = (o.xBand ? d3.scalePoint().domain(o.xBand).padding(0.3) : d3.scaleLinear().domain(o.xDomain || d3.extent(allx))).range([0, f.iw]);
    if (!o.xBand && !o.xDomain) x.nice();
    const y = d3.scaleLinear().domain(o.yDomain || pad(d3.extent(ally), 0.08)).range([f.ih, 0]).nice();
    axes(f, x, y, { xLabel: o.xLabel, yLabel: o.yLabel, grid: o.grid !== false, yFormat: o.yFormat, xTicks: o.xTicks, yLabelOffset: 22,
      xAxis: o.xBand ? d3.axisBottom(x).tickSize(2.5).tickPadding(2) : undefined });
    if (o.diagonal) f.g.append('line').attr('x1', x(0)).attr('y1', y(0)).attr('x2', x(1)).attr('y2', y(1)).attr('stroke', RULE).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    if (o.hline != null) f.g.append('line').attr('x1', 0).attr('x2', f.iw).attr('y1', y(o.hline)).attr('y2', y(o.hline)).attr('stroke', RULE).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    const line = d3.line().x(p => x(p[0])).y(p => y(p[1])).curve(o.step ? d3.curveStepAfter : d3.curveLinear);
    o.series.forEach(s => {
      if (s.band) f.g.append('path').attr('d', d3.area().x(p => x(p[0])).y0(p => y(p[1])).y1(p => y(p[2]))(s.band)).attr('fill', s.color).attr('fill-opacity', 0.15);
      f.g.append('path').attr('d', line(s.points)).attr('fill', 'none').attr('stroke', s.color).attr('stroke-width', s.width || 1.1).attr('stroke-dasharray', s.dash || null);
      if (s.dots) f.g.selectAll(null).data(s.points).join('circle').attr('cx', p => x(p[0])).attr('cy', p => y(p[1])).attr('r', 2).attr('fill', s.color).attr('stroke', '#fff').attr('stroke-width', 0.4)
        .on('mousemove', (ev, p) => tip.show(`<b>${esc(s.label)}</b><br>${esc(o.xLabel || 'x')}: ${p[0]}<br>${esc(o.yLabel || 'y')}: ${fmt(p[1], 3)}`, ev)).on('mouseleave', tip.hide);
    });
    if (o.legend !== false) legend(f, o.series.map(s => ({ label: s.label, color: s.color, type: 'line', dash: s.dash })), { y: -8 });
    // invisible hover rects for crosshair on band axes
    return f;
  }

  // ---------- grouped bar with error bars ----------
  function barsH(el, o) { // horizontal bars, items [{label, value, err, color}]
    const n = o.items.length, rowH = o.rowH || 12;
    const f = frame(el, { width: o.width, height: 22 + n * rowH + 24, margin: { t: 10, l: o.labelWidth || 80, b: 26, r: 22 } });
    const x = d3.scaleLinear().domain(o.domain || [0, d3.max(o.items, d => d.value + (d.err || 0)) * 1.1]).range([0, f.iw]).nice();
    const y = d3.scaleBand().domain(o.items.map(d => d.label)).range([0, f.ih]).paddingInner(0.28);
    const xa = f.g.append('g').attr('transform', `translate(0,${f.ih})`).call(d3.axisBottom(x).ticks(o.ticks || 5, o.xFormat).tickSize(2.5).tickPadding(2).tickSizeOuter(0)); styleAxis(xa);
    f.g.append('line').attr('x1', x(Math.max(x.domain()[0], Math.min(0, x.domain()[1])))).attr('x2', x(Math.max(x.domain()[0], Math.min(0, x.domain()[1])))).attr('y1', 0).attr('y2', f.ih).attr('stroke', INK).attr('stroke-width', LW.axis);
    if (o.vline != null) f.g.append('line').attr('x1', x(o.vline)).attr('x2', x(o.vline)).attr('y1', 0).attr('y2', f.ih).attr('stroke', MUTED).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    f.g.append('text').attr('x', f.iw / 2).attr('y', f.ih + 20).attr('text-anchor', 'middle').attr('font-size', FS.label).text(o.xLabel);
    const rows = f.g.selectAll('g.r').data(o.items).join('g').attr('class', 'r').attr('transform', d => `translate(0,${y(d.label)})`)
      .on('mousemove', (ev, d) => tip.show(o.tooltip ? o.tooltip(d) : `<b>${esc(d.label)}</b><br>${fmt(d.value, 3)}`, ev)).on('mouseleave', tip.hide);
    const x0 = x(Math.max(x.domain()[0], Math.min(0, x.domain()[1])));
    rows.append('rect').attr('x', d => Math.min(x0, x(d.value))).attr('height', y.bandwidth()).attr('width', d => Math.abs(x(d.value) - x0)).attr('fill', d => d.color).attr('rx', 0.8);
    rows.filter(d => d.err).append('line').attr('x1', d => x(d.value - d.err)).attr('x2', d => x(d.value + d.err)).attr('y1', y.bandwidth() / 2).attr('y2', y.bandwidth() / 2).attr('stroke', INK).attr('stroke-width', 0.6);
    rows.append('text').attr('x', -4).attr('y', y.bandwidth() / 2 + 2.1).attr('text-anchor', 'end').attr('font-size', FS.tick).text(d => d.label);
    rows.append('text').attr('x', d => (d.value < 0 ? x0 + 3 : x(d.value + (d.err || 0)) + 3)).attr('y', y.bandwidth() / 2 + 2.1).attr('font-size', FS.note).attr('fill', INK2).text(d => (o.valueFormat ? o.valueFormat(d.value) : fmt(d.value, 2)));
    return f;
  }

  // ---------- enrichment dot plot ----------
  function dotplot(el, o) {
    const n = o.items.length, rowH = 10.5;
    const f = frame(el, { width: o.width || 'onehalf', height: 30 + n * rowH + 28, margin: { t: 14, l: o.labelWidth || 150, b: 28, r: 62 } });
    const x = d3.scaleLinear().domain(o.xDomain || [0, d3.max(o.items, d => d.x) * 1.12]).range([0, f.iw]).nice();
    const y = d3.scaleBand().domain(o.items.map(d => d.id)).range([0, f.ih]).padding(0.25);
    const smax = d3.max(o.items, d => d.size) || 1, rS = d3.scaleSqrt().domain([0, smax]).range([0, Math.min(5, y.step() * 0.48)]);
    const cmax = Math.max(2, d3.max(o.items, d => d.color) || 2), col = sequential([0, cmax], 'blue');
    axes(f, x, y, { xLabel: o.xLabel, yAxis: d3.axisLeft(y).tickSize(0).tickPadding(4).tickFormat(id => { const it = o.items.find(d => d.id === id); const s = it ? it.name : id; return s.length > (o.maxChars || 46) ? s.slice(0, (o.maxChars || 46) - 1) + '…' : s; }), grid: false, xTicks: 4 });
    f.g.selectAll('.tick line').filter(function () { return false; });
    f.g.append('g').selectAll('line').data(o.items).join('line').attr('x1', 0).attr('x2', f.iw).attr('y1', d => y(d.id) + y.bandwidth() / 2).attr('y2', d => y(d.id) + y.bandwidth() / 2).attr('stroke', GRID).attr('stroke-width', 0.4);
    if (o.vline != null) f.g.append('line').attr('x1', x(o.vline)).attr('x2', x(o.vline)).attr('y1', 0).attr('y2', f.ih).attr('stroke', MUTED).attr('stroke-width', LW.thin).attr('stroke-dasharray', '2,1.5');
    f.g.selectAll('circle.d').data(o.items).join('circle').attr('class', 'd').attr('cx', d => x(d.x)).attr('cy', d => y(d.id) + y.bandwidth() / 2).attr('r', d => Math.max(1.2, rS(d.size)))
      .attr('fill', d => col(d.color)).attr('stroke', INK2).attr('stroke-width', 0.35)
      .on('mousemove', (ev, d) => tip.show(o.tooltip(d), ev)).on('mouseleave', tip.hide);
    // legends: size + colour
    const lg = f.g.append('g').attr('transform', `translate(${f.iw + 12},0)`);
    lg.append('text').attr('y', 0).attr('font-size', FS.note).text(o.sizeTitle || 'Genes');
    const svals = [...new Set([Math.max(1, Math.round(smax / 4)), Math.max(1, Math.round(smax / 2)), smax])];
    svals.forEach((v, i) => { lg.append('circle').attr('cx', 5).attr('cy', 9 + i * 10).attr('r', rS(v)).attr('fill', 'none').attr('stroke', INK2).attr('stroke-width', 0.4); lg.append('text').attr('x', 13).attr('y', 11 + i * 10).attr('font-size', FS.note).text(v); });
    const gy = 9 + svals.length * 10 + 10, id = 'dg' + Math.random().toString(36).slice(2, 8);
    const grad = f.svg.append('defs').append('linearGradient').attr('id', id).attr('x1', 0).attr('x2', 0).attr('y1', 1).attr('y2', 0);
    d3.range(0, 1.0001, 0.1).forEach(t => grad.append('stop').attr('offset', t).attr('stop-color', col(t * cmax)));
    lg.append('text').attr('y', gy).attr('font-size', FS.note).text(o.colorTitle || '−log₁₀ q');
    lg.append('rect').attr('y', gy + 4).attr('width', 5).attr('height', 34).attr('fill', `url(#${id})`);
    const ls = d3.scaleLinear().domain([0, cmax]).range([gy + 38, gy + 4]);
    styleAxis(lg.append('g').attr('transform', 'translate(5,0)').call(d3.axisRight(ls).ticks(3).tickSize(2).tickPadding(1.5))).select('.domain').remove();
    return f;
  }

  // ---------- network (deterministic force layout) ----------
  function network(el, o) {
    const f = frame(el, { width: o.width || 'onehalf', aspect: o.aspect || 0.7, margin: { t: 14, l: 6, b: 6, r: 6 } });
    const nodes = o.nodes.map(d => ({ ...d })), links = o.links.map(d => ({ ...d }));
    // seed positions on a circle per view to make layout deterministic
    const viewsSet = [...new Set(nodes.map(nd => nd.view))];
    nodes.forEach((nd, i) => { const a = (viewsSet.indexOf(nd.view) / Math.max(1, viewsSet.length)) * 2 * Math.PI + (i % 7) * 0.25; const r = nd.type === 'factor' ? 0 : f.ih * 0.35; nd.x = f.iw / 2 + r * Math.cos(a) + (i % 3); nd.y = f.ih / 2 + r * Math.sin(a) + (i % 5); });
    const sim = d3.forceSimulation(nodes).randomSource(d3.randomLcg(0.557))
      .force('link', d3.forceLink(links).id(d => d.id).distance(l => (l.kind === 'factor' ? 48 : 34)).strength(0.6))
      .force('charge', d3.forceManyBody().strength(-48)).force('collide', d3.forceCollide(d => (d.type === 'factor' ? 12 : 7)))
      .force('x', d3.forceX(f.iw / 2).strength(0.05)).force('y', d3.forceY(f.ih / 2).strength(0.07)).stop();
    for (let i = 0; i < 400; i++) sim.tick();
    // fit into frame
    const xs = d3.extent(nodes, d => d.x), ys = d3.extent(nodes, d => d.y), padd = 22;
    const sx = (f.iw - 2 * padd) / Math.max(1, xs[1] - xs[0]), sy = (f.ih - 2 * padd) / Math.max(1, ys[1] - ys[0]), s = Math.min(sx, sy, 1.6);
    const ox = f.iw / 2 - s * (xs[0] + xs[1]) / 2, oy = f.ih / 2 - s * (ys[0] + ys[1]) / 2;
    nodes.forEach(nd => { nd.px = ox + s * nd.x; nd.py = oy + s * nd.y; });
    const wmax = d3.max(links, l => Math.abs(l.value)) || 1;
    const lk = f.g.append('g').selectAll('line').data(links).join('line').attr('x1', l => l.source.px).attr('y1', l => l.source.py).attr('x2', l => l.target.px).attr('y2', l => l.target.py)
      .attr('stroke', l => (l.value >= 0 ? o.posColor : o.negColor)).attr('stroke-opacity', 0.7).attr('stroke-width', l => 0.35 + 1.6 * Math.abs(l.value) / wmax)
      .on('mousemove', (ev, l) => tip.show(o.linkTooltip(l), ev)).on('mouseleave', tip.hide);
    const nd = f.g.append('g').selectAll('g').data(nodes).join('g').attr('transform', d => `translate(${d.px},${d.py})`).style('cursor', 'pointer')
      .on('mousemove', (ev, d) => tip.show(o.nodeTooltip(d), ev)).on('mouseleave', tip.hide)
      .on('click', (ev, d) => { lk.attr('stroke-opacity', l => (l.source.id === d.id || l.target.id === d.id ? 0.95 : 0.12)); o.onClick && o.onClick(d); });
    nd.filter(d => d.type === 'factor').append('rect').attr('x', -10).attr('y', -6).attr('width', 20).attr('height', 12).attr('rx', 2).attr('fill', INK);
    nd.filter(d => d.type === 'factor').append('text').attr('text-anchor', 'middle').attr('y', 2.2).attr('font-size', 6).attr('font-weight', 600).attr('fill', '#fff').text(d => d.label);
    nd.filter(d => d.type !== 'factor').append('path').attr('d', d => symbolPath(d.shape || 'circle', 3.4)).attr('fill', d => d.color).attr('stroke', '#fff').attr('stroke-width', 0.6);
    nd.filter(d => d.type !== 'factor').append('text').attr('x', 4.6).attr('y', -3.2).attr('font-size', 5.8).attr('fill', INK).text(d => (d.label.length > 18 ? d.label.slice(0, 17) + '…' : d.label));
    legend(f, o.legendItems, { y: -8, align: 'left', x: 0 });
    return f;
  }

  // ---------- export ----------
  function serialize(svgEl) {
    const clone = svgEl.cloneNode(true);
    const wmm = svgEl.getAttribute('data-width-mm'), hmm = svgEl.getAttribute('data-height-mm');
    clone.setAttribute('width', wmm + 'mm'); clone.setAttribute('height', hmm + 'mm'); clone.removeAttribute('style');
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
  }
  function download(name, blob) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); }
  function exportSVG(svgEl, name) { download(name + '.svg', new Blob([serialize(svgEl)], { type: 'image/svg+xml' })); }
  function exportPNG(svgEl, name, dpi) {
    dpi = dpi || 600; const wmm = +svgEl.getAttribute('data-width-mm'), hmm = +svgEl.getAttribute('data-height-mm');
    const W = Math.round(wmm / 25.4 * dpi), H = Math.round(hmm / 25.4 * dpi);
    const img = new Image(); const url = URL.createObjectURL(new Blob([serialize(svgEl)], { type: 'image/svg+xml' }));
    img.onload = () => { const c = document.createElement('canvas'); c.width = W; c.height = H; const ctx = c.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H); ctx.drawImage(img, 0, 0, W, H); URL.revokeObjectURL(url); c.toBlob(b => download(`${name}_${dpi}dpi.png`, setDPI ? b : b), 'image/png'); };
    img.src = url;
  }
  const setDPI = true;

  root.Charts = { WIDTHS, MM, FONT, tip, fmt, fmtP, esc, scatter, volcano, sampleBars, sampleBoxes, histogram, lollipop, heatmap, strip, lines, barsH, dotplot, network, diverging, sequential, exportSVG, exportPNG, serialize };
})(window);
