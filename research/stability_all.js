/* Monodromy matrices for every orbit in a claims file, in the format stability.py reads.
 * Usage: node stability_all.js results/claims.json runs/claims_monodromy.json
 *        python3 stability.py runs/claims_monodromy.json > results/claims_stability.json */
"use strict";
const fs = require("fs");
const S = require("./slice.js");
const { monodromy } = require("./monodromy.js");
const [inp, out] = process.argv.slice(2);
if (!inp || !out) { console.error("Usage: node stability_all.js claims.json out_monodromy.json"); process.exit(1); }
const res = [];
for (const c of JSON.parse(fs.readFileSync(inp, "utf8"))) {
  const { m, x, v } = S.ic(c.m3, c.v1, c.v2);
  const y0 = [x[0], x[1], x[3], x[4], x[6], x[7], v[0], v[1], v[3], v[4], v[6], v[7]];
  res.push({ id: c.id, T: c.T, m, y0, M: monodromy(c.m3, c.v1, c.v2, c.T) });
  console.error(c.id);
}
fs.writeFileSync(out, JSON.stringify(res));
