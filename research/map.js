/* Return-distance map of the isosceles slice for the findings page: for every (v1, v2) on an
 * N x N grid, the smallest phase-space return distance min_t |Y(t) - Y(0)| over 0.5 < t < Tmax,
 * and the time at which it happens. Periodic orbits are the dark points (d -> 0).
 * Usage: node map.js m3 N Tmax out.json */
"use strict";
const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const os = require("os"), fs = require("fs"), path = require("path");
const { NBody, classify } = require(path.join(__dirname, "../web/solver.js"));
const S = require("./slice.js");

function scan(m3, v1, v2, Tmax) {
  const { m, x, v } = S.ic(m3, v1, v2);
  const y0 = [...x, ...v];
  const sim = new NBody(m, x, v, { tol: 1e-10 });
  let best = Infinity, tBest = 0, fate = 0; // fate 1: escaped (stopped early)
  const buf = new Float64Array(9), bv = new Float64Array(9);
  try {
    while (sim.t < Tmax && sim.steps < 40000) {
      sim.step();
      const a = sim.prev[18];
      for (let k = 1; k <= 4; k++) {
        const t = a + ((sim.t - a) * k) / 4;
        if (t <= 0.5 || t > Tmax) continue;
        sim.sample(t, buf, bv);
        let s = 0;
        for (let i = 0; i < 9; i++) { const p = buf[i] - y0[i], q = bv[i] - y0[9 + i]; s += p * p + q * q; }
        if (s < best) { best = s; tBest = t; }
      }
      if ((sim.steps & 15) === 0 && classify(sim.x, sim.v, m, 1, 6).state === "escape") { fate = 1; break; }
    }
  } catch (e) { fate = 2; }
  return [Math.sqrt(best), tBest, fate];
}

if (isMainThread) {
  const [m3s, Ns, Tmaxs, out] = process.argv.slice(2);
  const m3 = +m3s, N = +Ns, Tmax = +Tmaxs, vm = S.vmax(m3);
  const d = new Array(N * N).fill(-1), T = new Array(N * N).fill(0), F = new Array(N * N).fill(0);
  let next = 0, done = 0, alive = os.cpus().length; const t0 = Date.now();
  for (let w = 0; w < os.cpus().length; w++) {
    const worker = new Worker(__filename, { workerData: { m3, N, Tmax, vm } });
    const feed = () => worker.postMessage(next < N ? next++ : -1);
    worker.on("message", ({ r, row }) => {
      for (let c = 0; c < N; c++) { const [dd, tt, ff] = row[c]; d[r * N + c] = dd; T[r * N + c] = tt; F[r * N + c] = ff; }
      if (++done % Math.max(1, Math.floor(N / 10)) === 0) console.log(`${Math.round((100 * done) / N)}%  ${((Date.now() - t0) / 1000).toFixed(0)}s`);
      feed();
    });
    worker.on("exit", () => {
      if (--alive) return;
      fs.writeFileSync(out, JSON.stringify({ m3, N, Tmax, vmax: vm, d: d.map((q) => +q.toPrecision(3)), T: T.map((q) => +q.toFixed(2)), fate: F }));
      console.log(`done in ${((Date.now() - t0) / 1000).toFixed(0)}s`);
    });
    feed();
  }
} else {
  const { m3, N, Tmax, vm } = workerData;
  parentPort.on("message", (r) => {
    if (r < 0) return process.exit(0);
    const v2 = (vm * (r + 0.5)) / N, row = [];
    for (let c = 0; c < N; c++) {
      const v1 = (vm * (c + 0.5)) / N;
      row.push(S.energy(m3, v1, v2) >= 0 ? [-1, 0, 0] : scan(m3, v1, v2, Tmax));
    }
    parentPort.postMessage({ r, row });
  });
}
