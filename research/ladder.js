/* Extend the ladder A(Ba)^k B a(bA)^k b at fixed m3 by extrapolating (v1, v2, T) from the previous
 * rungs and refining. A rung counts only if Newton converges AND its word is exactly the predicted
 * one. Usage: node ladder.js m3 kmax out.json seeds.jsonl
 *   seeds.jsonl: a process.js *_orbits.jsonl holding the lower rungs found by the grid search.
 *   Also writes local box scans to runs/ladder_box_<m3>_k<k>.jsonl. */
"use strict";
const fs = require("fs"), path = require("path");
const { NBody } = require(path.join(__dirname, "../web/solver.js"));
const S = require("./slice.js"), TP = require("./topology.js"), R = require("./refine.js");
const CONV = { pole: 1, flip: false, up: 1 };
const word = (k) => "A" + "Ba".repeat(k) + "B" + "a" + "bA".repeat(k) + "b";
function wordOf(m3, p) {
  const { m, x, v } = S.ic(m3, p[0], p[1]);
  return TP.canonical(TP.wordFrom(TP.crossings(NBody, m, x, v, p[2], { tol: 1e-13 }), CONV) || "");
}
/* the predicted rung lies outside Newton's basin: scan a small box around the prediction for
   near-returns at about the predicted period, then refine the best few */
function localScan(m3, guess, k, prev) {
  for (const [N, shrink] of [[120, 1], [240, 0.5], [320, 0.25]]) {
    const f = localScanOnce(m3, guess, k, prev, N, shrink);
    if (f) return f;
  }
  return null;
}
function localScanOnce(m3, guess, k, prev, N, shrink) {
  const { execFileSync } = require("child_process");
  const dv1 = shrink * Math.max(1.5e-3, Math.abs(guess[0] - prev[0])), dv2 = shrink * Math.max(5e-4, Math.abs(guess[1] - prev[1]));
  const tmp = path.join(__dirname, "runs", `ladder_box_${m3}_k${k}.jsonl`);
  execFileSync("node", [path.join(__dirname, "search.js"), String(m3), String(N), String(Math.ceil(guess[2] * 1.12)), tmp, "0.25",
    String(guess[0] - dv1), String(guess[0] + dv1), String(guess[1] - dv2), String(guess[1] + dv2)], { stdio: "ignore" });
  const hits = [];
  for (const l of fs.readFileSync(tmp, "utf8").trim().split("\n").slice(1)) {
    const c = JSON.parse(l);
    for (const [t, d] of c.hits) if (Math.abs(t - guess[2]) < 0.1 * guess[2]) hits.push([d, c.v1, c.v2, t]);
  }
  hits.sort((x, y) => x[0] - y[0]);
  for (const h of hits.slice(0, 8)) {
    const f = R.refine(m3, [h[1], h[2], h[3]], { maxIt: 40 });
    if (f.ok && wordOf(m3, f.p) === TP.canonical(word(k))) return f;
  }
  return null;
}
const [m3s, kmaxs, out, seedFile] = process.argv.slice(2);
if (!seedFile) { console.error("Usage: node ladder.js m3 kmax out.json seeds.jsonl (a process.js *_orbits.jsonl)"); process.exit(1); }
const m3 = +m3s, kmax = +kmaxs;
// seeds: orbits from a grid-search result file whose word matches the ladder
const seeds = fs.readFileSync(seedFile, "utf8").trim().split("\n").map(JSON.parse)
  .map((o) => { for (let k = 1; k <= 30; k++) if (o.canon === TP.canonical(word(k))) return { k, p: [o.v1, o.v2, o.T] }; return null; })
  .filter(Boolean).sort((a, b) => a.k - b.k);
// several orbits can share a low-k word (distinct branches); keep the branch whose period grows
// steadily with k: start from the longest-k seed and walk down, keeping the seed whose period is
// closest to the arithmetic prediction
const byK = new Map(); for (const s of seeds) (byK.get(s.k) || byK.set(s.k, []).get(s.k)).push(s.p);
const ks = [...byK.keys()].sort((a, b) => b - a);
const rungs = new Map([[ks[0], byK.get(ks[0])[0]]]);
for (let i = 1; i < ks.length; i++) {
  const k = ks[i], hi = rungs.get(k + 1), hi2 = rungs.get(k + 2);
  if (!hi) break;
  const want = hi2 ? 2 * hi[2] - hi2[2] : null;
  const pick = byK.get(k).filter((p) => (want ? Math.abs(p[2] - want) < 0.15 * want : p[2] < hi[2]))
    .sort((a, b) => (want ? Math.abs(a[2] - want) - Math.abs(b[2] - want) : b[2] - a[2]))[0];
  if (!pick) break;
  rungs.set(k, pick);
}
console.log("seed rungs:", [...rungs.keys()].join(","));
for (let k = Math.max(...rungs.keys()) + 1; k <= kmax; k++) {
  const a = rungs.get(k - 1), b = rungs.get(k - 2), c = rungs.get(k - 3);
  if (!a || !b) break;
  // quadratic extrapolation in k when three rungs are known, linear otherwise
  const guess = a.map((q, i) => (c ? 3 * q - 3 * b[i] + c[i] : 2 * q - b[i]));
  let f = R.refine(m3, guess, { maxIt: 40 });
  if (!f.ok) f = localScan(m3, guess, k, a);
  if (!f || !f.ok) { console.log(`k=${k}: no convergence`); break; }
  const w = wordOf(m3, f.p), ok = w === TP.canonical(word(k));
  console.log(`k=${k}: v=(${f.p[0].toFixed(10)}, ${f.p[1].toFixed(10)}) T=${f.p[2].toFixed(6)} dT=${(f.p[2] - a[2]).toFixed(5)} res=${f.res.toExponential(1)} word ${ok ? "as predicted" : "DIFFERENT: " + w}`);
  if (!ok) break;
  rungs.set(k, f.p);
}
const rows = [...rungs.entries()].sort((p, q) => p[0] - q[0]).map(([k, p]) => ({ k, m3, v1: p[0], v2: p[1], T: p[2], E: S.energy(m3, p[0], p[1]), Tstar: p[2] * Math.pow(-S.energy(m3, p[0], p[1]), 1.5), Lf: 4 * k + 4, word: word(k) }));
fs.writeFileSync(out, JSON.stringify(rows, null, 1));
