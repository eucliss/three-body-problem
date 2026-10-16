/* Refine scan candidates into exact periodic orbits, classify them topologically and compare with
 * the published catalogues.
 * Usage: node process.js candidates.jsonl out.jsonl [dmax] [maxPerPoint]
 * Each output line: {m3, v1, v2, T, res, E, Tstar, word, canon, Lf, minSep, power, known, matches} */
"use strict";
const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const os = require("os"), fs = require("fs"), path = require("path");
const { NBody } = require(path.join(__dirname, "../web/solver.js"));
const S = require("./slice.js"), TP = require("./topology.js"), R = require("./refine.js");
const CONV = { pole: 1, flip: false, up: 1 }; // calibrated against 12 published words
const CONVS = [CONV, { pole: 1, flip: true, up: 1 }, { pole: 1, flip: false, up: -1 }, { pole: 1, flip: true, up: -1 }];

/* primitive root of a cyclic word: w = u^k → (u, k) */
function primitive(w) {
  for (let L = 1; L <= w.length / 2; L++) if (w.length % L === 0 && w.slice(0, L).repeat(w.length / L) === w) return [w.slice(0, L), w.length / L];
  return [w, 1];
}

function classifyOrbit(m3, v1, v2, T) {
  const { m, x, v } = S.ic(m3, v1, v2);
  const cr = TP.crossings(NBody, m, x, v, T, { tol: 1e-13 });
  const words = CONVS.map((c) => TP.wordFrom(cr, c));
  const w = words[0];
  // minimum pair separation over one period
  const sim = new NBody(m, x, v, { tol: 1e-12 });
  let minSep = Infinity; const b = new Float64Array(9);
  while (sim.t < T) {
    sim.step();
    for (let k = 1; k <= 8; k++) {
      sim.sample(sim.prev[18] + ((sim.t - sim.prev[18]) * k) / 8, b);
      for (const [i, j] of [[0, 1], [1, 2], [0, 2]]) minSep = Math.min(minSep, Math.hypot(b[3 * i] - b[3 * j], b[3 * i + 1] - b[3 * j + 1]));
    }
  }
  return { cr: cr.length, words, w, minSep };
}

if (isMainThread) {
  const [inp, out, dmaxs, maxPPs] = process.argv.slice(2);
  const dmax = +(dmaxs || 0.25), maxPP = +(maxPPs || 2);
  const lines = fs.readFileSync(inp, "utf8").trim().split("\n");
  const meta = JSON.parse(lines[0]).meta;
  const jobs = [];
  for (const l of lines.slice(1)) {
    const c = JSON.parse(l);
    for (const [t, d] of c.hits.slice(0, maxPP)) if (d < dmax) jobs.push([c.v1, c.v2, t, d]);
  }
  jobs.sort((a, b) => a[3] - b[3]); // most promising first
  const cat = JSON.parse(fs.readFileSync(path.join(__dirname, "data/catalog_words.json"), "utf8")).entries.filter((e) => e.set !== "free-fall");
  console.log(`m3=${meta.m3}: ${jobs.length} refinement jobs`);
  const fd = fs.openSync(out, "w");
  let next = 0, done = 0, ok = 0; const t0 = Date.now();
  const nW = os.cpus().length; let alive = nW;
  const seen = new Map();
  for (let w = 0; w < nW; w++) {
    const worker = new Worker(__filename, { workerData: { m3: meta.m3 } });
    const feed = () => worker.postMessage(next < jobs.length ? jobs[next++] : null);
    worker.on("message", (r) => {
      done++;
      if (r) {
        ok++;
        const canons = r.words.filter(Boolean).map((x) => TP.canonical(x));
        const matches = cat.filter((e) => canons.includes(e.canon)).map((e) => `${e.set}:${e.label}`);
        r.known = matches.length > 0; r.matches = matches.slice(0, 6);
        const key = r.canon + "|" + r.Tstar.toFixed(3);
        if (!seen.has(key)) { seen.set(key, true); fs.writeSync(fd, JSON.stringify(r) + "\n"); }
      }
      if (done % 200 === 0) console.log(`${done}/${jobs.length}  converged ${ok}  distinct ${seen.size}  ${((Date.now() - t0) / 1000).toFixed(0)}s`);
      feed();
    });
    worker.on("exit", () => { if (--alive === 0) { fs.closeSync(fd); console.log(`done: ${done} jobs, ${ok} converged, ${seen.size} distinct orbits, ${((Date.now() - t0) / 1000).toFixed(0)}s`); } });
    feed();
  }
} else {
  const { m3 } = workerData;
  parentPort.on("message", (job) => {
    if (!job) return process.exit(0);
    const [v1, v2, t] = job;
    let res = null;
    try {
      const f = R.refine(m3, [v1, v2, t], { maxIt: 30, patience: +(process.env.PATIENCE || 16) });
      if (f.ok && f.p[2] > 0.5 && S.energy(m3, f.p[0], f.p[1]) < 0) {
        const [p1, p2, T] = f.p;
        const c = classifyOrbit(m3, p1, p2, T);
        if (c.w !== null) {
          const [root, k] = primitive(TP.canonical(c.w));
          const E = S.energy(m3, p1, p2);
          res = { m3, v1: p1, v2: p2, T, res: f.res, E, Tstar: T * Math.pow(-E, 1.5), word: c.w, canon: TP.canonical(c.w), root, power: k,
            Lf: c.w.length, crossings: c.cr, minSep: c.minSep, words: c.words };
        }
      }
    } catch (e) { res = null; }
    parentPort.postMessage(res);
  });
}
