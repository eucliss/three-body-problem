/* Grid scan of the isosceles-collinear slice for near-periodic orbits.
 * For each (v1, v2) on an N x N grid over [0, vmax]^2 with E < 0, integrate to Tmax and record
 * the local minima of the phase-space return distance |Y(t) - Y(0)| (positions and velocities)
 * for t > tmin. Escapes stop early. Writes candidates (d < dthr) as JSON lines.
 * Usage: node search.js m3 N Tmax out.jsonl [dthr] [v1lo v1hi v2lo v2hi]
 */
"use strict";
const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const os = require("os"), fs = require("fs"), path = require("path");
const { NBody, classify } = require(path.join(__dirname, "../web/solver.js"));
const S = require("./slice.js");

function scan(m3, v1, v2, Tmax, dthr) {
  const { m, x, v } = S.ic(m3, v1, v2);
  const y0 = [...x, ...v];
  const sim = new NBody(m, x, v, { tol: 1e-11 });
  const tmin = 0.5;
  const dist = (y) => { let s = 0; for (let i = 0; i < 18; i++) { const q = y[i] - y0[i]; s += q * q; } return Math.sqrt(s); };
  const hits = [];
  let dPrev2 = Infinity, dPrev = Infinity, tPrev = 0;
  const buf = new Float64Array(9), bv = new Float64Array(9), yy = new Float64Array(18);
  try {
    while (sim.t < Tmax && sim.steps < 60000) {
      sim.step();
      // sample the return distance a few times per step (Hermite), enough to catch minima
      const a = sim.prev[18];
      for (let k = 1; k <= 4; k++) {
        const t = a + ((sim.t - a) * k) / 4;
        sim.sample(t, buf, bv);
        for (let i = 0; i < 9; i++) { yy[i] = buf[i]; yy[9 + i] = bv[i]; }
        const d = dist(yy);
        if (t > tmin && dPrev < dPrev2 && dPrev <= d && dPrev < dthr) hits.push([+tPrev.toFixed(6), +dPrev.toPrecision(4)]);
        dPrev2 = dPrev; dPrev = d; tPrev = t;
      }
      if ((sim.steps & 15) === 0) {
        const c = classify(sim.x, sim.v, m, 1, 6);
        if (c.state === "escape") break;
      }
    }
  } catch (e) { /* stalled: keep what we have */ }
  hits.sort((p, q) => p[1] - q[1]);
  return hits.slice(0, 4);
}

if (isMainThread) {
  const [m3s, Ns, Tmaxs, out, dthrs, ...box] = process.argv.slice(2);
  const m3 = +m3s, N = +Ns, Tmax = +Tmaxs, dthr = +(dthrs || 0.25);
  const vm = S.vmax(m3);
  const [v1lo, v1hi, v2lo, v2hi] = box.length === 4 ? box.map(Number) : [0, vm, 0, vm];
  const rows = [...Array(N).keys()];
  const nW = os.cpus().length;
  const fd = fs.openSync(out, "w");
  let next = 0, done = 0, nCand = 0;
  const t0 = Date.now();
  fs.writeSync(fd, JSON.stringify({ meta: { m3, N, Tmax, dthr, box: [v1lo, v1hi, v2lo, v2hi] } }) + "\n");
  let alive = nW;
  for (let w = 0; w < nW; w++) {
    const worker = new Worker(__filename, { workerData: { m3, N, Tmax, dthr, box: [v1lo, v1hi, v2lo, v2hi] } });
    const feed = () => worker.postMessage(next < N ? rows[next++] : -1);
    worker.on("message", (msg) => {
      for (const c of msg) { fs.writeSync(fd, JSON.stringify(c) + "\n"); nCand++; }
      done++;
      if (done % Math.max(1, Math.floor(N / 20)) === 0) console.log(`${Math.round((100 * done) / N)}%  ${((Date.now() - t0) / 1000).toFixed(0)}s  candidates ${nCand}`);
      feed();
    });
    worker.on("exit", () => { if (--alive === 0) { fs.closeSync(fd); console.log(`done m3=${m3} N=${N} Tmax=${Tmax} in ${((Date.now() - t0) / 1000).toFixed(0)}s, ${nCand} candidate rows`); } });
    feed();
  }
} else {
  const { m3, N, Tmax, dthr, box } = workerData;
  parentPort.on("message", (r) => {
    if (r < 0) return process.exit(0);
    const v2 = box[2] + ((box[3] - box[2]) * (r + 0.5)) / N;
    const res = [];
    for (let c = 0; c < N; c++) {
      const v1 = box[0] + ((box[1] - box[0]) * (c + 0.5)) / N;
      if (S.energy(m3, v1, v2) >= 0) continue;
      const hits = scan(m3, v1, v2, Tmax, dthr);
      if (hits.length) res.push({ v1: +v1.toFixed(9), v2: +v2.toFixed(9), hits });
    }
    parentPort.postMessage(res);
  });
}
