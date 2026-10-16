/* Continuation of a periodic-orbit branch parametrized by its period T: unknowns (v1, v2, m3),
 * residual Y(T) - Y(0). Unlike continuation in m3, this passes smoothly through turning points
 * in m3 (folds), which show up as m3 reversing direction while T keeps growing.
 * Usage: node contT.js name m3 v1 v2 T dT Tstop out.json */
"use strict";
const fs = require("fs"), path = require("path");
const R = require("./refine.js"), S = require("./slice.js"), C = require("./continue.js");
const norm = (r) => Math.sqrt(r.reduce((s, q) => s + q * q, 0));
function resid(q, T, tol) { if (!(q[2] > 0.005) || S.energy(q[2], q[0], q[1]) >= 0) return null; return R.residual(q[2], [q[0], q[1], T], tol); }
function solve3(A, b) {
  const M = A.map((row, i) => [...row, b[i]]);
  for (let c = 0; c < 3; c++) { let p = c; for (let r = c + 1; r < 3; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r; [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < 3; r++) { const f = M[r][c] / M[c][c]; for (let k = c; k < 4; k++) M[r][k] -= f * M[c][k]; } }
  const x = [0, 0, 0]; for (let r = 2; r >= 0; r--) { let s = M[r][3]; for (let k = r + 1; k < 3; k++) s -= M[r][k] * x[k]; x[r] = s / M[r][r]; } return x;
}
/* Levenberg–Marquardt on (v1, v2, m3) at fixed T */
function refineAtT(q0, T) {
  let q = q0.slice(), tol = 1e-12, r = resid(q, T, tol); if (!r) return null;
  let f = norm(r), mu = 1e-3;
  for (let it = 0; it < 30; it++) {
    if (f < 1e-8 && tol > 1e-13) { tol = 1e-13; r = resid(q, T, tol); if (!r) return null; f = norm(r); }
    if (f < 1e-11 && tol <= 1e-13) return { q, res: f };
    const J = [];
    for (let k = 0; k < 3; k++) { const h = 1e-8, qq = q.slice(); qq[k] += h; const rk = resid(qq, T, tol); if (!rk) return null; J.push(rk.map((v, i) => (v - r[i]) / h)); }
    const JTJ = [0, 1, 2].map((a) => [0, 1, 2].map((b) => J[a].reduce((s, _, i) => s + J[a][i] * J[b][i], 0)));
    const JTr = [0, 1, 2].map((a) => J[a].reduce((s, _, i) => s + J[a][i] * r[i], 0));
    let improved = false;
    for (let tries = 0; tries < 12; tries++) {
      const dq = solve3(JTJ.map((row, a) => row.map((v, b) => (a === b ? v * (1 + mu) : v))), JTr.map((v) => -v));
      const qn = q.map((v, k) => v + dq[k]); const rn = resid(qn, T, tol);
      if (rn && norm(rn) < f) { q = qn; r = rn; f = norm(rn); mu = Math.max(mu / 10, 1e-12); improved = true; break; }
      mu *= 10;
    }
    if (!improved) break;
  }
  return f < 1e-10 ? { q, res: f } : null;
}
function contT(q0, T0, dT0, Tstop) {
  const pts = [{ T: T0, q: q0 }];
  let dT = dT0;
  while (Math.abs(dT) > 1e-7 && (dT0 > 0 ? pts[pts.length - 1].T < Tstop : pts[pts.length - 1].T > Tstop)) {
    const a = pts[pts.length - 1], b = pts[pts.length - 2];
    const Tn = a.T + dT;
    const guess = b ? a.q.map((v, i) => v + ((v - b.q[i]) * (Tn - a.T)) / (a.T - b.T)) : a.q;
    const s = refineAtT(guess, Tn);
    if (s && Math.abs(s.q[2] - a.q[2]) < 0.05) { pts.push({ T: Tn, q: s.q, res: s.res }); dT = Math.sign(dT0) * Math.min(Math.abs(dT) * 1.4, Math.abs(dT0)); }
    else dT /= 2;
  }
  return pts;
}
if (require.main === module) {
  const [name, m3, v1, v2, T, dT, Tstop, out] = process.argv.slice(2);
  const pts = contT([+v1, +v2, +m3], +T, +dT, +Tstop);
  const ms = pts.map((p) => p.q[2]);
  const turns = []; for (let i = 1; i < ms.length - 1; i++) if ((ms[i] - ms[i - 1]) * (ms[i + 1] - ms[i]) < 0) turns.push({ m3: +ms[i].toFixed(6), T: +pts[i].T.toFixed(4) });
  const last = pts[pts.length - 1];
  const w0 = C.wordAt(+m3, [+v1, +v2, +T]), w1 = C.wordAt(last.q[2], [last.q[0], last.q[1], last.T]);
  fs.writeFileSync(out, JSON.stringify(pts.map((p) => ({ T: p.T, v1: p.q[0], v2: p.q[1], m3: p.q[2] }))));
  console.log(JSON.stringify({ name, points: pts.length, m3range: [Math.min(...ms), Math.max(...ms)].map((q) => +q.toFixed(5)), turningPoints: turns, last: { T: +last.T.toFixed(4), m3: +last.q[2].toFixed(6) }, sameWord: w0 === w1 }));
}
module.exports = { contT, refineAtT };
