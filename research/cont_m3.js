/* Natural continuation in m3 (secant predictor in v1, v2, T; Newton at fixed m3), used to push a
 * branch past the point where T-continuation stalls. Usage: node cont_m3.js branch.json dm m3stop out.json */
"use strict";
const fs = require("fs"), R = require("./refine.js"), S = require("./slice.js"), TP = require("./topology.js");
const { NBody } = require("../web/solver.js");
const [file, dms, stops, out] = process.argv.slice(2);
const br = JSON.parse(fs.readFileSync(file, "utf8")).sort((a, b) => a.m3 - b.m3);
let a = br[br.length - 2], b = br[br.length - 1], dm = +dms;
const pts = [];
const word = (p) => { const { m, x, v } = S.ic(p.m3, p.v1, p.v2); return TP.canonical(TP.wordFrom(TP.crossings(NBody, m, x, v, p.T, { tol: 1e-13 }), { pole: 1, flip: false, up: 1 })); };
const w0 = word(b);
while (b.m3 < +stops && Math.abs(dm) > 1e-5) {
  const m3 = b.m3 + dm, f = dm / (b.m3 - a.m3);
  const guess = [b.v1 + f * (b.v1 - a.v1), b.v2 + f * (b.v2 - a.v2), b.T + f * (b.T - a.T)];
  const r = R.refine(m3, guess, { maxIt: 25 });
  const ok = r.ok && Math.abs(r.p[2] - guess[2]) < 0.05 * guess[2] && Math.hypot(r.p[0] - guess[0], r.p[1] - guess[1]) < 0.01;
  if (!ok) { dm /= 2; continue; }
  const c = { m3, v1: r.p[0], v2: r.p[1], T: r.p[2], res: r.res };
  const w = word(c);
  if (w !== w0) { console.log("word changed at", m3, w); dm /= 2; continue; }
  pts.push(c); a = b; b = c;
  if (pts.length % 10 === 0) { console.log(JSON.stringify(c)); fs.writeFileSync(out, JSON.stringify(pts)); }
  dm = Math.min(dm * 1.3, +dms * 4);
}
fs.writeFileSync(out, JSON.stringify(pts));
console.log("stopped at", b.m3, "T", b.T, "points", pts.length);
