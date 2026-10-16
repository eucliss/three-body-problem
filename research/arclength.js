/* Pseudo-arclength continuation of a periodic orbit in (v1, v2, T, m3), which follows a branch
 * around folds (saddle-node bifurcations) where natural continuation in m3 stops.
 * Unknowns z = (v1, v2, T, m3); residual = [Y(T)-Y(0) (12), w * t·(z - z_pred)].
 * Scaled coordinates: (v1, v2, T/20, m3). Writes the branch to a JSON file.
 * Usage: node arclength.js name m3 v1 v2 T dirSign maxPoints ds out.json */
"use strict";
const fs = require("fs"), path = require("path");
const R = require("./refine.js"), S = require("./slice.js");
const SC = [1, 1, 1 / 20, 1];
const toS = (z) => z.map((q, i) => q * SC[i]), fromS = (s) => s.map((q, i) => q / SC[i]);
function res(z, tol) { if (!(z[3] > 0.005 && z[2] > 0.5)) return null; return R.residual(z[3], [z[0], z[1], z[2]], tol); }
const norm = (a) => Math.sqrt(a.reduce((s, q) => s + q * q, 0));

function solveLS(J, r) { // least squares J dx = -r via normal equations (4 unknowns)
  const n = J[0].length, A = [...Array(n)].map((_, a) => [...Array(n)].map((_, b) => J.reduce((s, row) => s + row[a] * row[b], 0)));
  const g = [...Array(n)].map((_, a) => -J.reduce((s, row, i) => s + row[a] * r[i], 0));
  for (let c = 0; c < n; c++) { let p = c; for (let q = c + 1; q < n; q++) if (Math.abs(A[q][c]) > Math.abs(A[p][c])) p = q; [A[c], A[p]] = [A[p], A[c]]; [g[c], g[p]] = [g[p], g[c]];
    for (let q = c + 1; q < n; q++) { const f = A[q][c] / A[c][c]; for (let k = c; k < n; k++) A[q][k] -= f * A[c][k]; g[q] -= f * g[c]; } }
  const x = Array(n).fill(0); for (let q = n - 1; q >= 0; q--) { let s = g[q]; for (let k = q + 1; k < n; k++) s -= A[q][k] * x[k]; x[q] = s / A[q][q]; }
  return x;
}

/* corrector: Gauss-Newton on 13 residuals in scaled coordinates */
function correct(zs0, tan, zsPred, tol) {
  let zs = zsPred.slice();
  for (let it = 0; it < 12; it++) {
    const z = fromS(zs), r0 = res(z, tol); if (!r0) return null;
    const ac = 10 * tan.reduce((s, t, i) => s + t * (zs[i] - zsPred[i]), 0);
    const r = [...r0, ac];
    if (norm(r0) < 1e-11) return { zs, res: norm(r0) };
    const J = r.map(() => []);
    for (let k = 0; k < 4; k++) {
      const h = 1e-8, zp = zs.slice(); zp[k] += h;
      const rp = res(fromS(zp), tol); if (!rp) return null;
      const acp = 10 * tan.reduce((s, t, i) => s + t * (zp[i] - zsPred[i]), 0);
      [...rp, acp].forEach((q, i) => { J[i][k] = (q - r[i]) / h; });
    }
    const dx = solveLS(J, r);
    zs = zs.map((q, i) => q + dx[i]);
    if (norm(dx) > 0.2) return null;
  }
  const r0 = res(fromS(zs), tol);
  return r0 && norm(r0) < 1e-10 ? { zs, res: norm(r0) } : null;
}

function trace(z0, dirSign, maxPoints, ds0, z1) {
  // initial tangent: a given second point, or a small natural-continuation step in m3
  if (!z1) {
    const f1 = R.refine(z0[3] + dirSign * 0.002, [z0[0], z0[1], z0[2]], { maxIt: 25 });
    if (!f1.ok) throw new Error("could not take the first step");
    z1 = [...f1.p, z0[3] + dirSign * 0.002];
  }
  let pts = [toS(z0), toS(z1)];
  let ds = ds0;
  const out = pts.map((p) => fromS(p));
  while (out.length < maxPoints && ds > 1e-5) {
    const a = pts[pts.length - 2], b = pts[pts.length - 1];
    let tan = b.map((q, i) => q - a[i]); const tn = norm(tan); tan = tan.map((q) => q / tn);
    const pred = b.map((q, i) => q + ds * tan[i]);
    const c = correct(b, tan, pred, 1e-12);
    if (c) {
      pts.push(c.zs); out.push(fromS(c.zs)); ds = Math.min(ds * 1.3, ds0);
    } else ds /= 2;
    const z = out[out.length - 1];
    if (z[3] < 0.01 || z[3] > 3 || z[2] > 400) break;
  }
  return out.map((z) => ({ v1: z[0], v2: z[1], T: z[2], m3: z[3], E: S.energy(z[3], z[0], z[1]) }));
}

if (require.main === module) {
  const [name, m3, v1, v2, T, dir, maxP, ds, out, seedBranch] = process.argv.slice(2);
  let z0 = [+v1, +v2, +T, +m3], z1 = null;
  if (seedBranch) { // start from two consecutive points of a natural-continuation branch
    const b = JSON.parse(fs.readFileSync(seedBranch, "utf8")).branch;
    z0 = [b[0].v1, b[0].v2, b[0].T, b[0].m3]; z1 = [b[1].v1, b[1].v2, b[1].T, b[1].m3];
  }
  const br = trace(z0, +dir, +maxP, +ds, z1);
  fs.writeFileSync(out, JSON.stringify({ name, branch: br }));
  const ms = br.map((p) => p.m3), turns = [];
  for (let i = 1; i < ms.length - 1; i++) if ((ms[i] - ms[i - 1]) * (ms[i + 1] - ms[i]) < 0) turns.push({ m3: +ms[i].toFixed(6), T: +br[i].T.toFixed(4) });
  console.log(JSON.stringify({ name, points: br.length, m3range: [Math.min(...ms), Math.max(...ms)].map((q) => +q.toFixed(5)), folds: turns, end: br[br.length - 1] }));
}
module.exports = { trace };
