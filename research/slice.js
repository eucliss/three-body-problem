/* The isosceles-collinear initial-condition slice used by Šuvakov & Dmitrašinović (2013) and the
 * Li & Liao catalogues: m1 = m2 = 1, r1 = (-1,0), r2 = (1,0), r3 = (0,0),
 * v1 = v2 = (v1, v2), v3 = -2 (v1, v2) / m3. Zero momentum and zero angular momentum. */
"use strict";
function ic(m3, v1, v2) {
  return {
    m: [1, 1, m3],
    x: [-1, 0, 0, 1, 0, 0, 0, 0, 0],
    v: [v1, v2, 0, v1, v2, 0, (-2 * v1) / m3, (-2 * v2) / m3, 0],
  };
}
/* total energy for the slice: K = |v|^2 (1 + 2/m3), U = 1/2 + 2 m3 */
const energy = (m3, v1, v2) => (v1 * v1 + v2 * v2) * (1 + 2 / m3) - (0.5 + 2 * m3);
const vmax = (m3) => Math.sqrt((0.5 + 2 * m3) / (1 + 2 / m3));
module.exports = { ic, energy, vmax };
