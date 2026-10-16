/* Exact orbit of a T-continuation branch at a target mass: secant root-finding on the period T
 * (each evaluation is a fixed-T solve that stays on the branch) until m3(T) = target. */
"use strict";
const { refineAtT } = require("./contT.js");
const S = require("./slice.js"), C = require("./continue.js");
function atMass(branch, target) {
  const out = [];
  for (let i = 1; i < branch.length; i++) {
    const a = branch[i - 1], b = branch[i];
    if ((a.m3 - target) * (b.m3 - target) > 0) continue;
    let T0 = a.T, T1 = b.T, q0 = [a.v1, a.v2, a.m3], q1 = [b.v1, b.v2, b.m3];
    for (let it = 0; it < 40 && Math.abs(q1[2] - target) > 1e-12; it++) {
      const Tn = T1 + ((target - q1[2]) * (T1 - T0)) / (q1[2] - q0[2]);
      const s = (Tn - T1) / (T1 - T0), guess = q1.map((v, k) => v + s * (v - q0[k]));
      const r = refineAtT(guess, Tn);
      if (!r) break;
      T0 = T1; q0 = q1; T1 = Tn; q1 = r.q;
    }
    if (Math.abs(q1[2] - target) < 1e-9) {
      const E = S.energy(q1[2], q1[0], q1[1]);
      out.push({ m3: target, v1: q1[0], v2: q1[1], T: T1, E, Tstar: T1 * Math.pow(-E, 1.5), word: C.wordAt(q1[2], [q1[0], q1[1], T1]) });
    }
  }
  return out;
}
module.exports = { atMass };
if (require.main === module) {
  const [file, ...targets] = process.argv.slice(2);
  const br = JSON.parse(require("fs").readFileSync(file, "utf8"));
  for (const t of targets.map(Number)) for (const o of atMass(br, t)) console.log(JSON.stringify(o));
}
