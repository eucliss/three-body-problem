/* Natural-parameter continuation of a periodic orbit in m3 (isosceles slice), with step halving to
 * locate where the branch ends (a fold). Each point is refined to residual < 1e-11 and its word is
 * re-computed, so the branch is verified to stay in one topological family.
 * Usage: node continue.js m3 v1 v2 T direction(+1|-1) [m3stop] [step] */
"use strict";
const path = require("path");
const { NBody } = require(path.join(__dirname, "../web/solver.js"));
const S = require("./slice.js"), TP = require("./topology.js"), R = require("./refine.js");
const CONV = { pole: 1, flip: false, up: 1 };

function wordAt(m3, p) {
  const { m, x, v } = S.ic(m3, p[0], p[1]);
  return TP.canonical(TP.wordFrom(TP.crossings(NBody, m, x, v, p[2], { tol: 1e-13 }), CONV) || "");
}

function continueBranch(m3, p, dir, m3stop, step0) {
  const canon0 = wordAt(m3, p);
  const branch = [{ m3, v1: p[0], v2: p[1], T: p[2], E: S.energy(m3, p[0], p[1]) }];
  let step = step0, prev = null, cur = { m3, p };
  while (step > 2e-6 && (dir > 0 ? cur.m3 < m3stop : cur.m3 > m3stop)) {
    const m3n = cur.m3 + dir * step;
    // secant predictor from the last two points
    let guess = cur.p.slice();
    if (prev) { const s = (m3n - cur.m3) / (cur.m3 - prev.m3); guess = cur.p.map((q, i) => q + s * (q - prev.p[i])); }
    const f = R.refine(m3n, guess, { maxIt: 25 });
    let ok = f.ok && S.energy(m3n, f.p[0], f.p[1]) < 0 && Math.abs(f.p[2] - guess[2]) < 0.2 * guess[2];
    if (ok && branch.length % 5 === 0) ok = wordAt(m3n, f.p) === canon0; // spot-check topology
    if (ok) {
      prev = cur; cur = { m3: m3n, p: f.p };
      branch.push({ m3: m3n, v1: f.p[0], v2: f.p[1], T: f.p[2], E: S.energy(m3n, f.p[0], f.p[1]), res: f.res });
      step = Math.min(step * 1.5, step0);
    } else step /= 2;
  }
  const last = branch[branch.length - 1];
  return { canon: canon0, endsAt: last.m3, reachedStop: dir > 0 ? last.m3 >= m3stop - 1e-9 : last.m3 <= m3stop + 1e-9, finalWord: wordAt(last.m3, [last.v1, last.v2, last.T]), branch };
}
module.exports = { continueBranch, wordAt };

if (require.main === module) {
  const [m3, v1, v2, T, dir, stop, st] = process.argv.slice(2).map(Number);
  const r = continueBranch(m3, [v1, v2, T], dir, stop ?? (dir > 0 ? 1 : 0.01), st ?? 0.01);
  console.log(JSON.stringify({ canon: r.canon, endsAt: r.endsAt, reachedStop: r.reachedStop, sameWordAtEnd: r.finalWord === r.canon, points: r.branch.length, last: r.branch[r.branch.length - 1] }));
}
