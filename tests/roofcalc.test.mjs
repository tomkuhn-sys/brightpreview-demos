// Run: node tests/roofcalc.test.mjs
// Builds synthetic houses + a bush on a fake DSM and checks the measurements.
import * as R from "../public/roofline/roofcalc.js";

let fails = 0;
const ok = (cond, msg) => { console.error((cond ? "  ok   " : "  FAIL ") + msg); if (!cond) fails++; };
const near = (a, b, tol, msg) => ok(Math.abs(a - b) <= tol, `${msg}: got ${a.toFixed(3)}, want ${b.toFixed(3)} ±${tol}`);

// ---- UTM
console.error("UTM");
{
  const z = R.utmZone(-122.0363); ok(z === 10, "Sunnyvale is zone 10");
  const p = R.toUTM(37.3688, -122.0363, z), q = R.fromUTM(p.e, p.n, z);
  near(q.lat, 37.3688, 1e-8, "round-trip lat"); near(q.lng, -122.0363, 1e-8, "round-trip lng");
  // 100 m east along the parallel ~ 0.0011305 deg lng at this latitude (geodesic); UTM grid scale ~0.99973 here
  const dLng = 100 / (Math.PI / 180 * 6378137 * Math.cos(37.3688 * Math.PI / 180) / Math.sqrt(1 - 0.00669438 * Math.sin(37.3688 * Math.PI / 180) ** 2));
  const p2 = R.toUTM(37.3688, -122.0363 + dLng, z);
  near(Math.hypot(p2.e - p.e, p2.n - p.n), 100 * 0.99973, 0.02, "100 m ground distance in grid meters");
  const s = R.toUTM(-33.9, 147.3, 55), s2 = R.fromUTM(s.e, s.n, 55, false);
  near(s2.lat, -33.9, 1e-8, "southern hemisphere round-trip");
  ok(Math.abs(R.convergence(37.3688, -122.0363, 10) - 0.58) < 0.05, "grid convergence ~0.58°");
}

// ---- synthetic scene
const px = 0.1, W = 600, H = 600, x0 = 585000, y0 = 4137060; // 60 m square
const dsm = { w: W, h: H, px, x0, y0, data: new Float32Array(W * H) };
const roof = new Uint8Array(W * H);
const rgb = [new Float32Array(W * H), new Float32Array(W * H), new Float32Array(W * H)];
const pitch = Math.atan(6 / 12) * 180 / Math.PI, tp = Math.tan(pitch * Math.PI / 180);
const groundZ = (e, n) => 20 + 0.01 * (e - x0); // gently sloping lot
// Hip-roof house: e in [x0+20, x0+32], n in [y0-40, y0-32] (12 x 8 m), eaves at 3 m
const hip = { e0: x0 + 20, e1: x0 + 32, n0: y0 - 40, n1: y0 - 32 };
// Gable garage: e in [x0+36, x0+46], n in [y0-40, y0-32], ridge E-W
const gab = { e0: x0 + 36, e1: x0 + 46, n0: y0 - 40, n1: y0 - 32 };
const bush = { e: x0 + 26, n: y0 - 44, r: 0.75, h: 1.2 };
const eave = 23; // absolute eave height (ground ~20 + 3)
for (let r = 0; r < H; r++) for (let c = 0; c < W; c++) {
  const { e, n } = R.px2utm(dsm, c, r), i = r * W + c;
  let z = groundZ(e, n);
  rgb[0][i] = 120; rgb[1][i] = 110; rgb[2][i] = 100; // grey ground
  if (e > hip.e0 && e < hip.e1 && n > hip.n0 && n < hip.n1) {
    const d = Math.min(e - hip.e0, hip.e1 - e, n - hip.n0, hip.n1 - n);
    z = eave + tp * d; roof[i] = 1; rgb[0][i] = 90; rgb[1][i] = 80; rgb[2][i] = 80;
  } else if (e > gab.e0 && e < gab.e1 && n > gab.n0 && n < gab.n1) {
    const d = Math.min(n - gab.n0, gab.n1 - n);
    z = eave + tp * d; roof[i] = 1; rgb[0][i] = 90; rgb[1][i] = 80; rgb[2][i] = 80;
  } else {
    const dd = Math.hypot(e - bush.e, n - bush.n);
    if (dd < bush.r) { z += bush.h * Math.sqrt(1 - (dd / bush.r) ** 2); rgb[0][i] = 50; rgb[1][i] = 120; rgb[2][i] = 40; }
  }
  dsm.data[i] = z + (Math.sin(i * 12.9898) * 43758.5453 % 1) * 0.04; // ±4 cm noise
}
// Solar-style segments (azimuth = compass direction the plane faces/downslope)
function seg(az, poly) {
  let ce = 0, cn = 0; for (const p of poly) { ce += p.e; cn += p.n; } ce /= poly.length; cn /= poly.length;
  const es = poly.map(p => p.e), ns = poly.map(p => p.n);
  // true plane height at centroid
  const facesN = az === 0, facesS = az === 180, facesE = az === 90;
  const box = poly === poly ? null : null;
  return { az, pitch, ce, cn, bbox: { e0: Math.min(...es), e1: Math.max(...es), n0: Math.min(...ns), n1: Math.max(...ns) }, area: 1, _h: null };
}
const P = (e, n) => ({ e, n });
const segs = [
  seg(180, [P(hip.e0, hip.n0), P(hip.e1, hip.n0), P(hip.e1 - 4, hip.n0 + 4), P(hip.e0 + 4, hip.n0 + 4)]), // south
  seg(0, [P(hip.e0, hip.n1), P(hip.e1, hip.n1), P(hip.e1 - 4, hip.n1 - 4), P(hip.e0 + 4, hip.n1 - 4)]),   // north
  seg(90, [P(hip.e1, hip.n0), P(hip.e1, hip.n1), P(hip.e1 - 4, hip.n1 - 4)]),                              // east
  seg(270, [P(hip.e0, hip.n0), P(hip.e0, hip.n1), P(hip.e0 + 4, hip.n1 - 4)]),                             // west
  seg(180, [P(gab.e0, gab.n0), P(gab.e1, gab.n0), P(gab.e1, gab.n0 + 4), P(gab.e0, gab.n0 + 4)]),          // garage south
  seg(0, [P(gab.e0, gab.n1), P(gab.e1, gab.n1), P(gab.e1, gab.n1 - 4), P(gab.e0, gab.n1 - 4)]),            // garage north
];
// plane height at center, deliberately in a different datum (+31.5 m) to exercise alignPlanes
const trueZ = (s) => {
  const inHip = s.ce < hip.e1 + 0.01;
  const b = inHip ? hip : gab;
  const d = inHip ? Math.min(s.ce - b.e0, b.e1 - s.ce, s.cn - b.n0, b.n1 - s.cn) : Math.min(s.cn - b.n0, b.n1 - s.cn);
  return eave + tp * d;
};
for (const s of segs) s.h0 = trueZ(s) + 31.5;

console.error("Roof planes");
const off = R.alignPlanes(dsm, segs);
near(off, -31.5, 0.05, "datum offset recovered");
const lab = R.labelFacets(dsm, roof, segs);
let labelled = 0, roofPx = 0; for (let i = 0; i < lab.length; i++) { if (roof[i]) roofPx++; if (lab[i] > 0) labelled++; }
ok(labelled / roofPx > 0.95, `pixels assigned to planes: ${(100 * labelled / roofPx).toFixed(1)}%`);

const { facets, edges } = R.buildEdges(dsm, lab, segs);
ok(facets.length === 6, `6 facets traced (got ${facets.length})`);
const front = { e: 0, n: -1 }; // street is to the south
R.assignSides(edges, { e: (hip.e0 + gab.e1) / 2, n: (hip.n0 + hip.n1) / 2 }, front);
const sum = (f) => edges.filter(f).reduce((a, e) => a + e.len, 0);
const rake3d = 4 / Math.cos(pitch * Math.PI / 180);
const hip3d = Math.hypot(4 * Math.SQRT2, 4 * tp);
console.table(edges.map(e => ({ type: e.type, side: e.side, facet: e.facet, len: +e.len.toFixed(2), slope: +e.slope.toFixed(1) })));
near(sum(e => e.type === "eave" && e.facet <= 3), 2 * 12 + 2 * 8, 0.6, "hip house eaves total (m)");
near(sum(e => e.type === "hip"), 4 * hip3d, 0.8, "hip lengths total (m)");
near(sum(e => e.type === "ridge" && e.facet <= 3), 4, 0.5, "hip house ridge (m)");
near(sum(e => e.type === "rake"), 4 * rake3d, 0.6, "garage rakes total, sloped (m)");
near(sum(e => e.type === "eave" && e.facet >= 4), 20, 0.5, "garage eaves (m)");
const rk = edges.filter(e => e.type === "rake");
ok(rk.every(e => Math.abs(e.slope - pitch) < 3), "rake slope ≈ roof pitch (26.6°)");
ok(edges.filter(e => e.type === "eave" && e.side === "front").length === 2, "two front eaves (house + garage)");
ok(rk.filter(e => e.side === "right").length === 2 && rk.filter(e => e.side === "left").length === 2, "garage rakes split left/right");

console.error("Presets");
R.presetSelect(edges, "front", { rakes: true });
near(sum(e => e.selected), 22, 0.6, "front preset = 12 + 10 m of eave");
R.presetSelect(edges, "frontSides", { rakes: true });
const fs = sum(e => e.selected);
near(fs, 22 + 8 + 8 + 4 * rake3d, 1.0, "front + sides (both hip-end eaves + all 4 garage rakes, which face left/right)");
R.presetSelect(edges, "all", { rakes: true, ridges: true, hips: true });
near(sum(e => e.selected), sum(() => true) - sum(e => e.type === "top" || e.type === "step"), 0.01, "all + ridges + hips selects everything except upper/step edges");

console.error("Front from roads");
{
  const c = { e: x0 + 30, n: y0 - 36 };
  const roads = [
    { name: "Maple Avenue", pts: [P(x0, y0 - 55), P(x0 + 60, y0 - 55)] },  // south, 19 m
    { name: "Oak Street", pts: [P(x0 + 52, y0), P(x0 + 52, y0 - 60)] },     // east, 22 m (corner lot)
  ];
  const r1 = R.frontFromRoads(c, roads, "Oak St");
  ok(r1 && r1.front.e > 0.99 && r1.matched, "corner lot: address street wins (faces east)");
  const r2 = R.frontFromRoads(c, roads, "");
  ok(r2 && r2.front.n < -0.99, "no street name: nearest road (faces south)");
}

console.error("Plants");
const plants = R.detectPlants(dsm, rgb, roof, { center: { e: x0 + 30, n: y0 - 36 } });
ok(plants.length === 1, `one plant found (got ${plants.length})`);
if (plants[0]) {
  near(plants[0].height, 1.2, 0.2, "bush height (m)");
  near(plants[0].diam, 1.5, 0.25, "bush width (m)");
  ok(plants[0].type === "bush" && plants[0].shape === "round", "classified as round bush");
}

console.error("Wrap math");
{
  const a = R.wrapLength({ shape: "round", height: 4, width: 5, spacingIn: 4 });
  near(a.feet, 149, 2, "4 ft tall x 5 ft wide round bush @4in");
  ok(a.rings === 12, "12 rings");
  const b = R.wrapLength({ shape: "round", height: 4, width: 5 / Math.PI, spacingIn: 4 });
  near(b.feet, 47, 2, "4 ft tall x 5 ft around");
  const c = R.wrapLength({ shape: "column", height: 6, width: 1, spacingIn: 6 });
  near(c.feet, 12 * Math.hypot(Math.PI, 0.5), 0.01, "column wrap = rings x spiral turn");
  const d = R.wrapLength({ shape: "box", height: 3, width: 2, length: 10, spacingIn: 6 });
  near(d.feet, 6 * Math.hypot(24, 0.5), 0.01, "hedge box wrap");
  const e = R.wrapLength({ shape: "cone", height: 6, width: 3, spacingIn: 6 });
  ok(e.feet > 0 && e.feet < R.wrapLength({ shape: "column", height: 6, width: 3, spacingIn: 6 }).feet * 0.55, "cone ~half of column");
}

console.error(fails ? `\n${fails} FAILED` : "\nall passed");
process.exit(fails ? 1 : 0);
