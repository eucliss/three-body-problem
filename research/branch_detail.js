/* Continue an orbit in m3 and save the full branch; report the orbit at a target mass. */
const fs = require("fs"), path = require("path");
const { NBody } = require(path.join(__dirname, "../web/solver.js"));
const C = require("./continue.js"), R = require("./refine.js"), S = require("./slice.js");
const [name, m3, v1, v2, T, dir, stop, target] = process.argv.slice(2);
const r = C.continueBranch(+m3, [+v1, +v2, +T], +dir, +stop, 0.01);
fs.writeFileSync(`runs/branch_${name}_${+dir > 0 ? "up" : "down"}.json`, JSON.stringify(r));
function minSep(m3, p) {
  const { m, x, v } = S.ic(m3, p[0], p[1]); const sim = new NBody(m, x, v, { tol: 1e-12 }); let ms = Infinity; const b = new Float64Array(9);
  while (sim.t < p[2]) { sim.step(); for (let k = 1; k <= 8; k++) { sim.sample(sim.prev[18] + (sim.t - sim.prev[18]) * k / 8, b); for (const [i, j] of [[0,1],[1,2],[0,2]]) ms = Math.min(ms, Math.hypot(b[3*i]-b[3*j], b[3*i+1]-b[3*j+1])); } }
  return ms;
}
const br = r.branch;
const near = br.reduce((a, b) => (Math.abs(b.m3 - +target) < Math.abs(a.m3 - +target) ? b : a));
const f = R.refine(+target, [near.v1, near.v2, near.T], { maxIt: 30 });
const last = br[br.length - 1], prev = br[Math.max(0, br.length - 6)];
console.log(JSON.stringify({ name, canon: r.canon, endsAt: last.m3,
  atTarget: f.ok ? { m3: +target, v1: f.p[0], v2: f.p[1], T: f.p[2], res: f.res, minSep: minSep(+target, f.p) } : "no convergence",
  end: { m3: last.m3, T: last.T, minSep: minSep(last.m3, [last.v1, last.v2, last.T]) },
  dT_dm3_near_end: (last.T - prev.T) / (last.m3 - prev.m3) }));
