/* Find the word convention that reproduces the published free-group words of the
   Šuvakov–Dmitrašinović orbits (Li & Liao 2017, Table a1), using our refined initial conditions. */
const path = require("path");
const { NBody } = require(path.join(__dirname, "../web/solver.js"));
const T = require("./topology.js");
const ORB = require(path.join(__dirname, "../data/periodic_orbits.json"));
const PUBLISHED = {
  figure_eight_SD: "abAB",
  butterfly_I: "(ab)2(AB)2",
  butterfly_II: "(ab)2(AB)2",
  bumblebee: "b2(ABab)2A2(baBA)2baB2(abAB)2a2(BAba)2BA",
  moth_I: "baBABabABA",
  butterfly_III: "(ab)2(ABA)(ba)2(BAB)",
  moth_II: "(abAB)2A(baBA)2B",
  moth_III: "ABabaBABabaBAbabABAbab",
  goggles: "(ab)2ABBA(ba)2BAAB",
  dragonfly: "b2(ABabAB)a2(BAbaBA)",
  yin_yang_Ia: "(ab)2ABAba(BAB)",
  yin_yang_Ib: "(ab)2ABAba(BAB)",
};
const results = {};
for (const [name, pub] of Object.entries(PUBLISHED)) {
  const o = ORB[name];
  const x = [-1, 0, 0, 1, 0, 0, 0, 0, 0], v = [o.p1, o.p2, 0, o.p1, o.p2, 0, -2 * o.p1, -2 * o.p2, 0];
  const cr = T.crossings(NBody, [1, 1, 1], x, v, o.T, { tol: 1e-13 });
  results[name] = { cr, pub: T.canonical(T.expand(pub)) };
}
for (const conv of T.CONVENTIONS) {
  let ok = 0; const bad = [];
  for (const [name, r] of Object.entries(results)) {
    const w = T.wordFrom(r.cr, conv);
    if (w !== null && T.canonical(w) === r.pub) ok++; else bad.push(`${name}:${w}`);
  }
  console.log(JSON.stringify(conv), `${ok}/${Object.keys(results).length}`, ok < 12 ? bad.slice(0, 3).join(" ") : "");
}
