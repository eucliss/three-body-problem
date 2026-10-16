/* Newton (Levenberg–Marquardt) refinement of a near-periodic orbit in the isosceles slice.
 * Unknowns (v1, v2, T); residual R = Y(T) - Y(0) over the planar positions and velocities.
 * Y(T) uses the solver's regularized exact-time landing, so T can be any real number. */
"use strict";
const path = require("path");
const { NBody } = require(path.join(__dirname, "../web/solver.js"));
const S = require("./slice.js");
const IDX = [0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16]; // planar x, y of positions then velocities

function flow(m3, v1, v2, T, tol) {
  const { m, x, v } = S.ic(m3, v1, v2);
  const sim = new NBody(m, x, v, { tol });
  let guard = 0;
  while (sim.t < T) { sim.step(); if (++guard > 400000) return null; }
  const out = new Float64Array(18);
  sim.stateAt(T, out);
  return { y: out, y0: [...x, ...v], minR: null };
}
function residual(m3, p, tol) {
  const f = flow(m3, p[0], p[1], p[2], tol);
  if (!f) return null;
  return IDX.map((i) => f.y[i] - f.y0[i]);
}
const norm = (r) => Math.sqrt(r.reduce((s, q) => s + q * q, 0));

function solve3(A, b) { // 3x3 linear solve (Gaussian elimination, partial pivoting)
  const M = A.map((row, i) => [...row, b[i]]);
  for (let c = 0; c < 3; c++) {
    let p = c; for (let r = c + 1; r < 3; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    if (Math.abs(M[c][c]) < 1e-300) return null;
    for (let r = c + 1; r < 3; r++) { const f = M[r][c] / M[c][c]; for (let k = c; k < 4; k++) M[r][k] -= f * M[c][k]; }
  }
  const xs = [0, 0, 0];
  for (let r = 2; r >= 0; r--) { let s = M[r][3]; for (let k = r + 1; k < 3; k++) s -= M[r][k] * xs[k]; xs[r] = s / M[r][r]; }
  return xs;
}

/* returns {p:[v1,v2,T], res, iters, ok} */
function refine(m3, p0, opts) {
  opts = opts || {};
  const maxIt = opts.maxIt ?? 40, target = opts.target ?? 1e-11;
  const patience = opts.patience ?? Infinity; // give up when not near-converged after this many iterations
  let p = p0.slice(), tol = 1e-11;
  let r = residual(m3, p, tol); if (!r) return { ok: false, p, res: Infinity, iters: 0 };
  let f = norm(r), mu = 1e-3;
  for (let it = 0; it < maxIt; it++) {
    if (f < 1e-8 && tol > 1e-13) { tol = 1e-13; r = residual(m3, p, tol); if (!r) break; f = norm(r); }
    if (f < target && tol <= 1e-13) return { ok: true, p, res: f, iters: it };
    if (it >= patience && f > 1e-6) break;
    // forward-difference Jacobian (12 x 3)
    const h = [1e-8, 1e-8, 1e-8 * Math.max(1, p[2])];
    const J = [];
    for (let k = 0; k < 3; k++) {
      const q = p.slice(); q[k] += h[k];
      const rk = residual(m3, q, tol); if (!rk) return { ok: false, p, res: f, iters: it };
      J.push(rk.map((v, i) => (v - r[i]) / h[k]));
    }
    // normal equations with LM damping
    const JTJ = [0, 1, 2].map((a) => [0, 1, 2].map((b) => J[a].reduce((s, _, i) => s + J[a][i] * J[b][i], 0)));
    const JTr = [0, 1, 2].map((a) => J[a].reduce((s, _, i) => s + J[a][i] * r[i], 0));
    let improved = false;
    for (let tries = 0; tries < 12; tries++) {
      const A = JTJ.map((row, a) => row.map((v, b) => (a === b ? v * (1 + mu) : v)));
      const dp = solve3(A, JTr.map((v) => -v));
      if (!dp) { mu *= 10; continue; }
      const q = p.map((v, k) => v + dp[k]);
      if (!(q[2] > 0.5)) { mu *= 10; continue; }
      const rq = residual(m3, q, tol);
      if (rq) {
        const fq = norm(rq);
        if (fq < f) { p = q; r = rq; f = fq; mu = Math.max(mu / 10, 1e-12); improved = true; break; }
      }
      mu *= 10;
    }
    if (!improved) break;
  }
  if (tol > 1e-13) { tol = 1e-13; const rr = residual(m3, p, tol); if (rr) f = norm(rr); }
  return { ok: f < target, p, res: f, iters: maxIt };
}
module.exports = { refine, residual, flow, IDX };
