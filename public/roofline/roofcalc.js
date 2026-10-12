// roofcalc.js — geometry for the Roofline & Plant Measure demo.
// Pure functions, no DOM, so they can be tested in Node.
// Coordinates: UTM meters (e = east, n = north). Grids are row-major, row 0 = north edge.

export const M2FT = 3.280839895;
const D2R = Math.PI / 180, R2D = 180 / Math.PI;

// ---------------------------------------------------------------- UTM (WGS84)
const A = 6378137, F = 1 / 298.257223563, E2 = F * (2 - F), EP2 = E2 / (1 - E2), K0 = 0.9996;
const cm = zone => ((zone - 1) * 6 - 180 + 3) * D2R;

export function utmZone(lng) { return Math.floor((lng + 180) / 6) + 1; }

export function toUTM(lat, lng, zone) {
  const phi = lat * D2R, lam = lng * D2R, l0 = cm(zone);
  const s = Math.sin(phi), c = Math.cos(phi), t = Math.tan(phi);
  const N = A / Math.sqrt(1 - E2 * s * s), T = t * t, C = EP2 * c * c, Aa = c * (lam - l0);
  const e4 = E2 * E2, e6 = e4 * E2;
  const M = A * ((1 - E2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * phi
    - (3 * E2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * Math.sin(2 * phi)
    + (15 * e4 / 256 + 45 * e6 / 1024) * Math.sin(4 * phi)
    - (35 * e6 / 3072) * Math.sin(6 * phi));
  const e = K0 * N * (Aa + (1 - T + C) * Aa ** 3 / 6 + (5 - 18 * T + T * T + 72 * C - 58 * EP2) * Aa ** 5 / 120) + 500000;
  let n = K0 * (M + N * t * (Aa * Aa / 2 + (5 - T + 9 * C + 4 * C * C) * Aa ** 4 / 24
    + (61 - 58 * T + T * T + 600 * C - 330 * EP2) * Aa ** 6 / 720));
  if (lat < 0) n += 10000000;
  return { e, n };
}

export function fromUTM(e, n, zone, north = true) {
  const x = e - 500000, y = north ? n : n - 10000000;
  const e4 = E2 * E2, e6 = e4 * E2;
  const mu = (y / K0) / (A * (1 - E2 / 4 - 3 * e4 / 64 - 5 * e6 / 256));
  const e1 = (1 - Math.sqrt(1 - E2)) / (1 + Math.sqrt(1 - E2));
  const p1 = mu + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * Math.sin(2 * mu) + (21 * e1 * e1 / 16 - 55 * e1 ** 4 / 32) * Math.sin(4 * mu)
    + (151 * e1 ** 3 / 96) * Math.sin(6 * mu) + (1097 * e1 ** 4 / 512) * Math.sin(8 * mu);
  const s = Math.sin(p1), c = Math.cos(p1), t = Math.tan(p1);
  const N1 = A / Math.sqrt(1 - E2 * s * s), T1 = t * t, C1 = EP2 * c * c, R1 = A * (1 - E2) / Math.pow(1 - E2 * s * s, 1.5);
  const D = x / (N1 * K0);
  const lat = p1 - (N1 * t / R1) * (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * EP2) * D ** 4 / 24
    + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * EP2 - 3 * C1 * C1) * D ** 6 / 720);
  const lng = cm(zone) + (D - (1 + 2 * T1 + C1) * D ** 3 / 6 + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * EP2 + 24 * T1 * T1) * D ** 5 / 120) / c;
  return { lat: lat * R2D, lng: lng * R2D };
}

// Angle (deg) from true north to grid north; subtract from a compass azimuth to get a grid azimuth.
export function convergence(lat, lng, zone) {
  return Math.atan(Math.tan(lng * D2R - cm(zone)) * Math.sin(lat * D2R)) * R2D;
}

// ---------------------------------------------------------------- grids
// grid: { w, h, data, x0, y0, px }  (x0,y0 = UTM of the top-left corner of pixel 0,0)
export function px2utm(g, col, row) { return { e: g.x0 + (col + 0.5) * g.px, n: g.y0 - (row + 0.5) * g.px }; }
export function utm2px(g, e, n) { return { col: (e - g.x0) / g.px, row: (g.y0 - n) / g.px }; }
export function cellAt(g, e, n) {
  const c = Math.floor((e - g.x0) / g.px), r = Math.floor((g.y0 - n) / g.px);
  return (c < 0 || r < 0 || c >= g.w || r >= g.h) ? -1 : r * g.w + c;
}

// Resample another grid (possibly other size/res) onto the base grid, nearest neighbour.
export function resampleTo(base, src, band = 0) {
  if (src.w === base.w && src.h === base.h && Math.abs(src.px - base.px) < 1e-9 &&
      Math.abs(src.x0 - base.x0) < 1e-6 && Math.abs(src.y0 - base.y0) < 1e-6) return src.bands ? src.bands[band] : src.data;
  const arr = src.bands ? src.bands[band] : src.data;
  const out = new Float32Array(base.w * base.h);
  for (let r = 0; r < base.h; r++) for (let c = 0; c < base.w; c++) {
    const p = px2utm(base, c, r), i = cellAt(src, p.e, p.n);
    out[r * base.w + c] = i < 0 ? NaN : arr[i];
  }
  return out;
}

// Separable sliding min/max over a square window of (2*rad+1) pixels. NaN treated as +Inf for min, -Inf for max.
export function slidingExtreme(data, w, h, rad, isMax) {
  const bad = isMax ? -Infinity : Infinity;
  const better = isMax ? (a, b) => a >= b : (a, b) => a <= b;
  const tmp = new Float32Array(w * h), out = new Float32Array(w * h);
  const line = (get, set, len) => {
    const dq = new Int32Array(len), vals = new Float64Array(len);
    for (let i = 0; i < len; i++) { const v = get(i); vals[i] = Number.isNaN(v) ? bad : v; }
    let head = 0, tail = 0, nextIn = 0;
    for (let i = 0; i < len; i++) {
      const hi = Math.min(len - 1, i + rad), lo = i - rad;
      while (nextIn <= hi) { while (tail > head && better(vals[nextIn], vals[dq[tail - 1]])) tail--; dq[tail++] = nextIn++; }
      while (dq[head] < lo) head++;
      set(i, vals[dq[head]]);
    }
  };
  for (let r = 0; r < h; r++) line(i => data[r * w + i], (i, v) => { tmp[r * w + i] = v; }, w);
  for (let c = 0; c < w; c++) line(i => tmp[i * w + c], (i, v) => { out[i * w + c] = v; }, h);
  return out;
}

function median(a) { if (!a.length) return 0; const s = [...a].sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
function percentile(a, p) { if (!a.length) return 0; const s = [...a].sort((x, y) => x - y); return s[Math.min(s.length - 1, Math.floor(p * (s.length - 1)))]; }

// ---------------------------------------------------------------- roof planes
// seg (from Solar API, already converted): { ce, cn, h0, pitch(deg), az(grid deg), bbox:{e0,n0,e1,n1}, area }
export function planeZ(seg, e, n) {
  const a = seg.az * D2R, slope = Math.tan(seg.pitch * D2R);
  return seg.h0 - slope * ((e - seg.ce) * Math.sin(a) + (n - seg.cn) * Math.cos(a));
}

// Shift plane heights so they line up with the DSM (handles any datum difference).
export function alignPlanes(dsm, segs) {
  const diffs = [];
  for (const s of segs) {
    const ring = [];
    for (let dy = -2; dy <= 2; dy++) for (let dx = -2; dx <= 2; dx++) {
      const i = cellAt(dsm, s.ce + dx * dsm.px, s.cn + dy * dsm.px);
      if (i >= 0 && !Number.isNaN(dsm.data[i])) ring.push(dsm.data[i]);
    }
    if (ring.length) diffs.push(median(ring) - s.h0);
  }
  const off = median(diffs);
  for (const s of segs) s.h0 += off;
  return off;
}

// Assign each roof pixel to the best-fitting roof plane. Returns Int32Array labels: -1 not roof, 0 roof but no plane, k+1 plane k.
export function labelFacets(dsm, roofMask, segs, opts = {}) {
  const tol = opts.tol ?? 0.45, pad = opts.pad ?? 0.8;
  const { w, h } = dsm;
  let lab = new Int32Array(w * h).fill(-1);
  const boxes = segs.map(s => ({ e0: s.bbox.e0 - pad, e1: s.bbox.e1 + pad, n0: s.bbox.n0 - pad, n1: s.bbox.n1 + pad }));
  for (let r = 0; r < h; r++) for (let c = 0; c < w; c++) {
    const i = r * w + c;
    if (!(roofMask[i] > 0)) continue;
    const z = dsm.data[i]; lab[i] = 0;
    if (Number.isNaN(z)) continue;
    const p = px2utm(dsm, c, r);
    let best = -1, bestErr = tol;
    for (let k = 0; k < segs.length; k++) {
      const b = boxes[k];
      if (p.e < b.e0 || p.e > b.e1 || p.n < b.n0 || p.n > b.n1) continue;
      const err = Math.abs(z - planeZ(segs[k], p.e, p.n));
      if (err < bestErr) { bestErr = err; best = k; }
    }
    if (best >= 0) lab[i] = best + 1;
  }
  for (let it = 0; it < 2; it++) lab = modeFilter(lab, w, h);
  keepLargestPerLabel(lab, w, h, segs.length);
  return lab;
}

function modeFilter(lab, w, h) {
  const out = new Int32Array(lab);
  const cnt = new Map();
  for (let r = 1; r < h - 1; r++) for (let c = 1; c < w - 1; c++) {
    const i = r * w + c; if (lab[i] < 0) continue;
    cnt.clear();
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      const v = lab[i + dy * w + dx]; if (v > 0) cnt.set(v, (cnt.get(v) || 0) + 1);
    }
    let bv = lab[i], bc = 0;
    for (const [v, n] of cnt) if (n > bc) { bc = n; bv = v; }
    if (bc >= 5) out[i] = bv;
  }
  return out;
}

export function components(isOn, w, h, conn8 = false) {
  const comp = new Int32Array(w * h).fill(-1), list = [];
  const stack = new Int32Array(w * h);
  for (let s = 0; s < w * h; s++) {
    if (comp[s] >= 0 || !isOn(s)) continue;
    const id = list.length, pix = []; let sp = 0; stack[sp++] = s; comp[s] = id;
    while (sp) {
      const i = stack[--sp]; pix.push(i);
      const r = (i / w) | 0, c = i - r * w;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        if (!dx && !dy) continue; if (!conn8 && dx && dy) continue;
        const rr = r + dy, cc = c + dx; if (rr < 0 || cc < 0 || rr >= h || cc >= w) continue;
        const j = rr * w + cc; if (comp[j] < 0 && isOn(j)) { comp[j] = id; stack[sp++] = j; }
      }
    }
    list.push(pix);
  }
  return { comp, list };
}

function keepLargestPerLabel(lab, w, h, nseg) {
  for (let k = 1; k <= nseg; k++) {
    const { list } = components(i => lab[i] === k, w, h);
    if (list.length < 2) continue;
    list.sort((a, b) => b.length - a.length);
    for (let j = 1; j < list.length; j++) for (const i of list[j]) lab[i] = 0;
  }
}

// ---------------------------------------------------------------- outlines
// Trace the outline of a pixel set along pixel boundaries; returns the largest loop as [{e,n}] in UTM.
export function traceOutline(isIn, g, bbox) {
  const { w } = g, W1 = w + 1;
  const [c0, r0, c1, r1] = bbox || [0, 0, g.w - 1, g.h - 1];
  const out = new Map(); // vertexKey -> array of [toKey]
  const add = (x0, y0, x1, y1) => { const k = x0 + y0 * W1; (out.get(k) || out.set(k, []).get(k)).push(x1 + y1 * W1); };
  const at = (c, r) => c >= 0 && r >= 0 && c < g.w && r < g.h && isIn(r * w + c);
  for (let r = r0; r <= r1; r++) for (let c = c0; c <= c1; c++) {
    if (!at(c, r)) continue;
    if (!at(c, r - 1)) add(c + 1, r, c, r);
    if (!at(c - 1, r)) add(c, r, c, r + 1);
    if (!at(c, r + 1)) add(c, r + 1, c + 1, r + 1);
    if (!at(c + 1, r)) add(c + 1, r + 1, c + 1, r);
  }
  let best = null, bestArea = 0;
  for (const [startKey, lst] of out) {
    while (lst.length) {
      const loop = []; let k = startKey;
      for (let guard = 0; guard < 1e6; guard++) {
        const nexts = out.get(k); if (!nexts || !nexts.length) break;
        const nk = nexts.pop(); loop.push(k); k = nk;
        if (k === startKey) break;
      }
      if (loop.length < 4) continue;
      const pts = loop.map(v => { const x = v % W1, y = (v / W1) | 0; return { e: g.x0 + x * g.px, n: g.y0 - y * g.px }; });
      const a = Math.abs(polyArea(pts));
      if (a > bestArea) { bestArea = a; best = pts; }
    }
  }
  return best;
}

export function polyArea(p) { let s = 0; for (let i = 0; i < p.length; i++) { const a = p[i], b = p[(i + 1) % p.length]; s += a.e * b.n - b.e * a.n; } return s / 2; }

function distToSeg(p, a, b) {
  const dx = b.e - a.e, dy = b.n - a.n, L = dx * dx + dy * dy;
  let t = L ? ((p.e - a.e) * dx + (p.n - a.n) * dy) / L : 0; t = Math.max(0, Math.min(1, t));
  return Math.hypot(p.e - (a.e + t * dx), p.n - (a.n + t * dy));
}
function dp(pts, eps) {
  if (pts.length < 3) return pts.slice();
  let idx = 0, dmax = 0;
  for (let i = 1; i < pts.length - 1; i++) { const d = distToSeg(pts[i], pts[0], pts[pts.length - 1]); if (d > dmax) { dmax = d; idx = i; } }
  if (dmax <= eps) return [pts[0], pts[pts.length - 1]];
  const a = dp(pts.slice(0, idx + 1), eps), b = dp(pts.slice(idx), eps);
  return a.slice(0, -1).concat(b);
}
export function simplifyClosed(pts, eps) {
  if (!pts || pts.length < 4) return pts;
  let far = 0, fd = 0;
  for (let i = 1; i < pts.length; i++) { const d = Math.hypot(pts[i].e - pts[0].e, pts[i].n - pts[0].n); if (d > fd) { fd = d; far = i; } }
  const a = dp(pts.slice(0, far + 1), eps), b = dp(pts.slice(far).concat([pts[0]]), eps);
  return a.slice(0, -1).concat(b.slice(0, -1));
}

// Snap edge directions to the building's main axes (multiples of 45°) and rebuild corners.
export function squareUp(poly, axisDeg, tolDeg = 12, minEdge = 0.3) {
  if (!poly || poly.length < 3) return poly;
  let edges = [];
  for (let i = 0; i < poly.length; i++) {
    const a = poly[i], b = poly[(i + 1) % poly.length];
    let ang = Math.atan2(b.n - a.n, b.e - a.e) * R2D;
    const rel = ((ang - axisDeg) % 45 + 45) % 45;
    const snap = rel < tolDeg ? -rel : (45 - rel < tolDeg ? 45 - rel : 0);
    ang += snap;
    edges.push({ m: { e: (a.e + b.e) / 2, n: (a.n + b.n) / 2 }, d: { e: Math.cos(ang * D2R), n: Math.sin(ang * D2R) }, a, b });
  }
  // merge consecutive near-parallel edges
  for (let changed = true; changed && edges.length > 3;) {
    changed = false;
    for (let i = 0; i < edges.length && edges.length > 3; i++) {
      const p = edges[i], q = edges[(i + 1) % edges.length];
      if (Math.abs(p.d.e * q.d.n - p.d.n * q.d.e) < 0.17 && p.d.e * q.d.e + p.d.n * q.d.n > 0) {
        const a = p.a, b = q.b, len = Math.hypot(b.e - a.e, b.n - a.n) || 1;
        const merged = { m: { e: (a.e + b.e) / 2, n: (a.n + b.n) / 2 }, d: { e: (b.e - a.e) / len, n: (b.n - a.n) / len }, a, b };
        if (i + 1 < edges.length) edges.splice(i, 2, merged);
        else { edges[i] = merged; edges.shift(); } // wrap-around: last edge + first edge
        changed = true; break;
      }
    }
  }
  const out = [];
  for (let i = 0; i < edges.length; i++) {
    const p = edges[(i - 1 + edges.length) % edges.length], q = edges[i];
    const cr = p.d.e * q.d.n - p.d.n * q.d.e;
    let v = q.a;
    if (Math.abs(cr) > 1e-6) {
      const t = ((q.m.e - p.m.e) * q.d.n - (q.m.n - p.m.n) * q.d.e) / cr;
      const x = { e: p.m.e + t * p.d.e, n: p.m.n + t * p.d.n };
      if (Math.hypot(x.e - q.a.e, x.n - q.a.n) < 2.0) v = x;
    }
    out.push(v);
  }
  // drop tiny edges
  for (let changed = true; changed && out.length > 3;) {
    changed = false;
    for (let i = 0; i < out.length; i++) {
      const a = out[i], b = out[(i + 1) % out.length];
      if (Math.hypot(b.e - a.e, b.n - a.n) < minEdge) { out.splice((i + 1) % out.length, 1); changed = true; break; }
    }
  }
  return out;
}

// Main axis of the building from the roof-plane azimuths (area-weighted circular mean mod 90°).
export function mainAxis(segs) {
  let sx = 0, sy = 0;
  for (const s of segs) { const a = s.az * 4 * D2R, w = s.area || 1; sx += w * Math.cos(a); sy += w * Math.sin(a); }
  const az = Math.atan2(sy, sx) * R2D / 4; // compass-style angle of downslope
  return 90 - az; // convert compass azimuth to math angle (from east, CCW)
}

// ---------------------------------------------------------------- edges
export const EDGE_TYPES = {
  eave: "Eave", rake: "Rake", ridge: "Ridge", hip: "Hip", valley: "Valley", top: "Upper edge", step: "Step"
};

// Returns { facets:[{k, poly}], edges:[...] }
export function buildEdges(dsm, lab, segs, opts = {}) {
  const px = dsm.px;
  const eps = opts.eps ?? Math.max(0.25, 1.5 * px);
  const axis = opts.axis ?? mainAxis(segs);
  const facets = [];
  for (let k = 0; k < segs.length; k++) {
    const s = segs[k];
    let n = 0; for (let i = 0; i < lab.length; i++) if (lab[i] === k + 1) n++;
    if (n * px * px < (opts.minFacet ?? 1.5)) continue;
    const raw = traceOutline(i => lab[i] === k + 1, dsm);
    if (!raw) continue;
    let poly = simplifyClosed(raw, eps);
    if (s.pitch > 2) poly = squareUp(poly, axis);
    if (poly && poly.length >= 3) facets.push({ k, poly, areaPlan: Math.abs(polyArea(poly)) });
  }
  const edges = [];
  const labAt = (e, n) => { const i = cellAt(dsm, e, n); return i < 0 ? -1 : lab[i]; };
  for (const f of facets) {
    const s = segs[f.k], P = f.poly, ccw = polyArea(P) > 0;
    const azr = s.az * D2R, dn = { e: Math.sin(azr), n: Math.cos(azr) }; // downslope dir
    for (let i = 0; i < P.length; i++) {
      const a = P[i], b = P[(i + 1) % P.length];
      const L = Math.hypot(b.e - a.e, b.n - a.n); if (L < 0.3) continue;
      const u = { e: (b.e - a.e) / L, n: (b.n - a.n) / L };
      const nrm = ccw ? { e: u.n, n: -u.e } : { e: -u.n, n: u.e }; // outward
      const tally = new Map(); let samples = 0;
      for (let t = 0.15; t <= 0.851; t += 0.1) {
        const q = { e: a.e + u.e * L * t + nrm.e * 0.45, n: a.n + u.n * L * t + nrm.n * 0.45 };
        const v = labAt(q.e, q.n); samples++;
        if (v > 0 && v !== f.k + 1) tally.set(v, (tally.get(v) || 0) + 1);
      }
      let other = 0, oc = 0; for (const [v, c] of tally) if (c > oc) { oc = c; other = v; }
      const internal = other > 0 && oc >= samples * 0.4;
      const j = internal ? other - 1 : -1;
      if (internal && j < f.k && facets.some(x => x.k === j)) continue; // the other facet reports it
      const along = Math.abs(u.e * dn.e + u.n * dn.n);   // 0 = level edge, 1 = straight up the slope
      const level = s.pitch < 5 || along < Math.sin(25 * D2R);
      const downhill = nrm.e * dn.e + nrm.n * dn.n;     // >0: edge is on the low side of the facet
      let type;
      if (!internal) type = s.pitch < 5 ? "eave" : level ? (downhill > 0 ? "eave" : "top") : "rake";
      else if (level) type = downhill < 0 ? "ridge" : "step";
      else type = downhill < 0 ? "hip" : "valley";
      const z = (p) => internal ? (planeZ(s, p.e, p.n) + planeZ(segs[j], p.e, p.n)) / 2 : planeZ(s, p.e, p.n);
      const za = z(a), zb = z(b), rise = Math.abs(zb - za);
      edges.push({
        id: edges.length, type, facet: f.k, other: j, a: { ...a, z: za }, b: { ...b, z: zb },
        normal: nrm, horiz: L, len: Math.hypot(L, rise), rise, slope: Math.atan2(rise, L) * R2D,
        pitch: s.pitch, mid: { e: (a.e + b.e) / 2, n: (a.n + b.n) / 2 }
      });
    }
  }
  return { facets, edges };
}

// Which side of the house an edge is on, given the unit "front" vector (house -> street).
export function sideOf(vec, front) {
  const L = Math.hypot(vec.e, vec.n) || 1, v = { e: vec.e / L, n: vec.n / L };
  const dot = v.e * front.e + v.n * front.n;
  if (dot >= Math.SQRT1_2) return "front";
  if (dot <= -Math.SQRT1_2) return "back";
  const left = { e: front.n, n: -front.e }; // viewer stands in the street facing the house
  return v.e * left.e + v.n * left.n > 0 ? "left" : "right";
}

export function assignSides(edges, center, front) {
  for (const ed of edges) {
    const outer = ed.other < 0 && ed.type !== "step";
    ed.side = sideOf(outer ? ed.normal : { e: ed.mid.e - center.e, n: ed.mid.n - center.n }, front);
  }
}

export const PRESET_SIDES = { front: ["front"], frontSides: ["front", "left", "right"], all: ["front", "left", "right", "back"], none: [] };
export function presetSelect(edges, preset, inc = {}) {
  const sides = PRESET_SIDES[preset] || [];
  for (const ed of edges) {
    const typeOk = ed.type === "eave" || (ed.type === "rake" && inc.rakes !== false) ||
      (ed.type === "ridge" && inc.ridges) || ((ed.type === "hip" || ed.type === "valley") && inc.hips);
    ed.selected = typeOk && sides.includes(ed.side);
  }
}

// ---------------------------------------------------------------- front of house from nearby roads
export function frontFromRoads(center, roads, streetName) {
  const norm = s => (s || "").toLowerCase().replace(/\b(street|st|avenue|ave|road|rd|drive|dr|court|ct|lane|ln|way|place|pl|boulevard|blvd|circle|cir|terrace|ter)\b\.?/g, "").replace(/[^a-z0-9]/g, "");
  const want = norm(streetName);
  let best = null;
  for (const rd of roads) {
    const named = want && norm(rd.name) && (norm(rd.name) === want);
    for (let i = 0; i < rd.pts.length - 1; i++) {
      const a = rd.pts[i], b = rd.pts[i + 1], dx = b.e - a.e, dy = b.n - a.n, L = dx * dx + dy * dy;
      let t = L ? ((center.e - a.e) * dx + (center.n - a.n) * dy) / L : 0; t = Math.max(0, Math.min(1, t));
      const p = { e: a.e + t * dx, n: a.n + t * dy }, d = Math.hypot(p.e - center.e, p.n - center.n);
      const score = d - (named ? 1000 : 0) + (rd.service ? 25 : 0);
      if (!best || score < best.score) best = { score, p, d, name: rd.name, matched: named };
    }
  }
  if (!best || best.d < 0.5) return null;
  const v = { e: best.p.e - center.e, n: best.p.n - center.n }, L = Math.hypot(v.e, v.n);
  return { front: { e: v.e / L, n: v.n / L }, road: best.name, matched: best.matched, dist: best.d };
}

// ---------------------------------------------------------------- plants
// dsm grid; rgb: [r,g,b] arrays on the dsm grid; roof: building mask on the dsm grid (any building).
export function detectPlants(dsm, rgb, roof, opts = {}) {
  const { w, h, px } = dsm;
  const center = opts.center, maxDist = opts.maxDist ?? 30;
  const ground = slidingExtreme(dsm.data, w, h, Math.round((opts.groundWin ?? 12) / 2 / px), false);
  const roofF = new Float32Array(w * h); for (let i = 0; i < w * h; i++) roofF[i] = roof[i] > 0 ? 1 : 0;
  const roofDil = slidingExtreme(roofF, w, h, Math.max(1, Math.round(0.6 / px)), true);
  const nd = new Float32Array(w * h);
  let veg = new Uint8Array(w * h);
  for (let i = 0; i < w * h; i++) {
    const z = dsm.data[i]; nd[i] = Number.isNaN(z) ? 0 : z - ground[i];
    if (roofDil[i] > 0) continue;
    const r = rgb[0][i], g = rgb[1][i], b = rgb[2][i], s = r + g + b + 1;
    const exg = (2 * g - r - b) / s;
    if ((nd[i] > (opts.minH ?? 0.3) && exg > 0.03) || (nd[i] > 2.5 && exg > -0.01)) veg[i] = 1;
  }
  // open (erode then dilate) to drop specks and thin fences
  const er = new Uint8Array(w * h);
  for (let r = 1; r < h - 1; r++) for (let c = 1; c < w - 1; c++) {
    const i = r * w + c; if (!veg[i]) continue;
    er[i] = veg[i - 1] & veg[i + 1] & veg[i - w] & veg[i + w];
  }
  veg = new Uint8Array(w * h);
  for (let r = 1; r < h - 1; r++) for (let c = 1; c < w - 1; c++) {
    const i = r * w + c; if (er[i] || er[i - 1] || er[i + 1] || er[i - w] || er[i + w]) veg[i] = 1;
  }
  const { comp, list } = components(i => veg[i] === 1, w, h, true);
  const plants = [];
  list.forEach((pix, id) => {
    const area = pix.length * px * px;
    if (area < (opts.minArea ?? 0.3)) return;
    let se = 0, sn = 0; const hs = [];
    for (const i of pix) { const r = (i / w) | 0, c = i - r * w, p = px2utm(dsm, c, r); se += p.e; sn += p.n; hs.push(nd[i]); }
    const ce = se / pix.length, cn = sn / pix.length;
    if (center && Math.hypot(ce - center.e, cn - center.n) > maxDist) return;
    let cee = 0, cnn = 0, cen = 0;
    for (const i of pix) { const r = (i / w) | 0, c = i - r * w, p = px2utm(dsm, c, r); const de = p.e - ce, dn = p.n - cn; cee += de * de; cnn += dn * dn; cen += de * dn; }
    cee /= pix.length; cnn /= pix.length; cen /= pix.length;
    const tr = cee + cnn, det = cee * cnn - cen * cen, disc = Math.sqrt(Math.max(0, tr * tr / 4 - det));
    const l1 = tr / 2 + disc, l2 = Math.max(tr / 2 - disc, 1e-6);
    const len = 4 * Math.sqrt(l1), wid = 4 * Math.sqrt(l2);
    const height = Math.max(0.3, percentile(hs, 0.9));
    const diam = 2 * Math.sqrt(area / Math.PI);
    let type = "bush", shape = "round";
    if (height >= 4 || area > 60) { type = "tree"; shape = "cone"; }
    else if (len / wid >= 2.5 && len >= 2.5) { type = "hedge"; shape = "box"; }
    let rs = h, re = 0, cs = w, ce2 = 0;
    for (const i of pix) { const r = (i / w) | 0, c = i - r * w; if (r < rs) rs = r; if (r > re) re = r; if (c < cs) cs = c; if (c > ce2) ce2 = c; }
    const outline = simplifyClosed(traceOutline(i => comp[i] === id, dsm, [cs, rs, ce2, re]), Math.max(0.15, px * 1.5));
    plants.push({
      e: ce, n: cn, area, height, diam, length: type === "hedge" ? len : diam, width: type === "hedge" ? wid : diam,
      type, shape, outline, dist: center ? Math.hypot(ce - center.e, cn - center.n) : 0
    });
  });
  plants.sort((a, b) => a.dist - b.dist);
  return plants.slice(0, opts.max ?? 40);
}

// ---------------------------------------------------------------- light wrapping math
// All inputs in feet except spacingIn (inches). Rings every `spacing` from bottom to top, each ring
// tilted up by one spacing (a spiral), length = sum of sqrt(perimeter^2 + spacing^2).
export function wrapLength({ shape, height, width, length, spacingIn = 4 }) {
  const s = spacingIn / 12;
  if (!(height > 0) || !(s > 0)) return { feet: 0, rings: 0 };
  const n = Math.max(1, Math.round(height / s)), dh = height / n;
  let total = 0;
  for (let k = 0; k < n; k++) {
    const y = (k + 0.5) * dh;
    let P;
    if (shape === "round") { const t = 2 * y / height - 1; P = Math.PI * width * Math.sqrt(Math.max(0, 1 - t * t)); }
    else if (shape === "cone") P = Math.PI * width * (1 - y / height);
    else if (shape === "box") P = 2 * ((length || width) + width);
    else P = Math.PI * width; // column
    total += Math.sqrt(P * P + dh * dh);
  }
  return { feet: total, rings: n };
}

export const toPitch12 = deg => 12 * Math.tan(deg * D2R);
