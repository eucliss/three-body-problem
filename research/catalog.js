/* Build an index of every published free-group word (canonical form, up to cyclic rotation and
 * inversion) from the Li & Liao catalogue markdown files, plus the Šuvakov–Dmitrašinović list.
 * Usage: node catalog.js <liao-repo-dir> out.json */
"use strict";
const fs = require("fs"), path = require("path");
const T = require("./topology.js");
const [dir, out] = process.argv.slice(2);
const strip = (s) => s.replace(/<[^>]*>/g, " ").replace(/&nbsp;/g, " ").replace(/\s+/g, " ").trim();
const entries = [];
for (const [file, set] of [["three-body-free-group-word.md", "equal-mass"], ["three-body-unequal-mass-free-group-word.md", "unequal-mass"], ["free-fall-3b-free-group-word.md", "free-fall"]]) {
  const p = path.join(dir, file);
  if (!fs.existsSync(p)) continue;
  const lines = fs.readFileSync(p, "utf8").split(/\n/).map(strip).filter(Boolean);
  let label = null;
  for (const l of lines) {
    if (/^[abAB]+$/.test(l)) { if (label) entries.push({ set, label, word: l, canon: T.canonical(l) }); label = null; }
    else if (/^(I|II|III|IV|V)\.[A-C]|^F\s*\d|i\.c\./.test(l)) label = l;
  }
}
const SD = { "figure-eight": "abAB", "butterfly I/II": "(ab)2(AB)2", bumblebee: "b2(ABab)2A2(baBA)2baB2(abAB)2a2(BAba)2BA", "moth I": "baBABabABA", "butterfly III": "(ab)2(ABA)(ba)2(BAB)", "moth II": "(abAB)2A(baBA)2B", "moth III": "ABabaBABabaBAbabABAbab", goggles: "(ab)2ABBA(ba)2BAAB", "butterfly IV": "((ba)2(AB)2)6A((ba)2(BA)2)6B", dragonfly: "b2(ABabAB)a2(BAbaBA)", "yin-yang I": "(ab)2ABAba(BAB)", "yin-yang II": "(abaBAB)3abaBAbab(ABAbab)3(AB)2" };
for (const [name, w] of Object.entries(SD)) entries.push({ set: "suvakov-dmitrasinovic-2013", label: name, word: T.expand(w), canon: T.canonical(T.expand(w)) });
const counts = {};
for (const e of entries) counts[e.set] = (counts[e.set] || 0) + 1;
fs.writeFileSync(out, JSON.stringify({ counts, entries }));
console.log(counts, "distinct canonical words:", new Set(entries.map((e) => e.canon)).size);
