/* Cluster scan candidates before refinement: keep the best few near-returns per (v1, v2, T) cell.
 * Usage: node prefilter.js in.jsonl out.jsonl [dmax] [cell] [tbin] [keep] [maxJobs]
 * Neighbouring grid points that return at the same time are almost always the same orbit, and
 * dense patches of near-misses near collision otherwise swamp the refinement queue. */
"use strict";
const fs = require("fs");
const [inp, out, dmaxs, cells, tbins, keeps, caps] = process.argv.slice(2);
const dmax = +(dmaxs || 0.1), cell = +(cells || 0.004), tbin = +(tbins || 0.5), keep = +(keeps || 3), cap = +(caps || Infinity);
const lines = fs.readFileSync(inp, "utf8").trim().split("\n");
const best = new Map();
for (const l of lines.slice(1)) {
  const c = JSON.parse(l);
  for (const [t, d] of c.hits) {
    if (d >= dmax) continue;
    const key = `${Math.floor(c.v1 / cell)},${Math.floor(c.v2 / cell)},${Math.round(t / tbin)}`;
    if (!best.has(key)) best.set(key, []);
    const b = best.get(key);
    b.push({ v1: c.v1, v2: c.v2, t, d }); b.sort((p, q) => p.d - q.d); if (b.length > keep) b.pop();
  }
}
const fd = fs.openSync(out, "w");
fs.writeSync(fd, lines[0] + "\n");
// rank: every cell's best candidate first (by d), then the cells' runners-up; truncate to the job budget
const ranked = [];
for (let r = 0; r < keep; r++) ranked.push(...[...best.values()].filter((bs) => bs[r]).map((bs) => bs[r]).sort((p, q) => p.d - q.d));
const kept = ranked.slice(0, cap);
let n = 0;
for (const b of kept) { n++; fs.writeSync(fd, JSON.stringify({ v1: b.v1, v2: b.v2, hits: [[b.t, b.d]] }) + "\n"); }
fs.closeSync(fd);
console.log(`${lines.length - 1} rows -> ${n} candidates in ${best.size} cells (dmax ${dmax}, cell ${cell}, tbin ${tbin}, keep ${keep})`);
