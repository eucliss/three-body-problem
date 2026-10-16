/* Parse the numeric tables (v1, v2, T, T*, Lf) of the Li–Liao "pictures" catalogue files.
 * Usage: node catalog_numeric.js <liao-repo-dir> out.json
 * Catalogue source: https://github.com/sjtu-liao/three-body at commit b96b019 (2026-02-15). */
"use strict";
const fs = require("fs"), path = require("path");
const [dir, out] = process.argv.slice(2);
const strip = (s) => s.replace(/<[^>]*>/g, " ").replace(/&nbsp;/g, " ").replace(/\s+/g, " ").trim();
const rows = [];
for (const [file, set, defaultM3] of [["three-body-pictures.md", "equal-mass", 1], ["three-body-unequal-mass-pictures.md", "unequal-mass", null]]) {
  const lines = fs.readFileSync(path.join(dir, file), "utf8").split(/\n/).map(strip).filter(Boolean);
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(/^((?:I|II|III|IV|V)\.[A-C])\s*i\.c\.\s*(\d+)\s*(?:\(([\d.]+)\))?$/);
    if (!m) continue;
    const nums = lines.slice(i + 1, i + 6).map(Number);
    if (nums.some((q) => !isFinite(q))) continue;
    rows.push({ set, label: `${m[1]} i.c. ${m[2]}`, m3: m[3] ? +m[3] : defaultM3, v1: nums[0], v2: nums[1], T: nums[2], Tstar: nums[3], Lf: nums[4] });
  }
}
fs.writeFileSync(out, JSON.stringify(rows));
const by = {}; for (const r of rows) by[`${r.set} m3=${r.m3}`] = (by[`${r.set} m3=${r.m3}`] || 0) + 1;
console.log(rows.length, by);
