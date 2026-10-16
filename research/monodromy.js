/* Monodromy matrix of a periodic orbit (planar, 12-dimensional phase space) by central
 * differences: column j = [Y(T; Y0 + h e_j) - Y(T; Y0 - h e_j)] / 2h, with exact-time landing.
 * Prints JSON {M: 12x12}. Eigenvalues are computed in Python (stability.py). */
"use strict";
const path = require("path");
const { NBody } = require(path.join(__dirname, "../web/solver.js"));
const S = require("./slice.js");
const IDX = [0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16];
function endState(m, y0, T) {
  const sim = new NBody(m, y0.slice(0, 9), y0.slice(9, 18), { tol: 1e-13 });
  while (sim.t < T) sim.step();
  const out = new Float64Array(18); sim.stateAt(T, out); return out;
}
function monodromy(m3, v1, v2, T, h) {
  h = h ?? 1e-7;
  const { m, x, v } = S.ic(m3, v1, v2);
  const y0 = [...x, ...v];
  const M = IDX.map(() => new Array(12).fill(0));
  IDX.forEach((jj, j) => {
    const yp = y0.slice(), ym = y0.slice(); yp[jj] += h; ym[jj] -= h;
    const ep = endState(m, yp, T), em = endState(m, ym, T);
    IDX.forEach((ii, i) => { M[i][j] = (ep[ii] - em[ii]) / (2 * h); });
  });
  return M;
}
module.exports = { monodromy };
if (require.main === module) {
  const [m3, v1, v2, T] = process.argv.slice(2).map(Number);
  const { m, x, v } = S.ic(m3, v1, v2);
  const y0 = [x[0], x[1], x[3], x[4], x[6], x[7], v[0], v[1], v[3], v[4], v[6], v[7]];
  console.log(JSON.stringify({ m3, v1, v2, T, m, y0, M: monodromy(m3, v1, v2, T) }));
}
