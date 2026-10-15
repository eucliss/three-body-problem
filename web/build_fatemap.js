/* Fate map of the equal-mass free-fall three-body problem.
 * Bodies 1, 2 start at (-1/2, 0), (1/2, 0); body 3 at (x, y); all at rest; G = m = 1.
 * Every pixel is a full integration with the same solver the page uses (web/solver.js) until
 * one body escapes or t = TMAX. Output: fatemap.json with one byte per pixel:
 *   255 = no escape by TMAX, else escaper * 64 + q, where q ∈ [0, 63] is log-scaled escape time.
 * Both halves are integrated directly (no mirroring), so every pixel can be replayed bit-for-bit
 * in the page by loading the same initial condition with the same tolerance.
 * Usage: node build_fatemap.js [W] [H] [TMAX]
 */
const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const os = require("os");
const fs = require("fs");
const path = require("path");
const { NBody, classify } = require("./solver.js");

const X_MAX = 1.0, Y_MAX = 1.2, TOL = 1e-11;
const T_LO = 2;

/* pixel centres; the page uses the identical expressions */
const pixelX = (c, W) => -X_MAX + (2 * X_MAX * (c + 0.5)) / W;
const pixelY = (r, H) => Y_MAX * (1 - (r + 0.5) / H);

function quantize(t, tmax) {
  const q = Math.round((63 * Math.log(Math.max(t, T_LO) / T_LO)) / Math.log(tmax / T_LO));
  return Math.max(0, Math.min(63, q));
}

/* Integrate one pixel. k = escaper (0-2), -1 = no escape by tmax, -2 = unresolved: the run hit
   the step cap, threw, or lost accuracy (energy error above 1e-6 over the 64 steps after the
   stop; checking a single instant would misfire at a pericentre, where T and U are huge and E is
   only known to ~U*tol). */
function fate(px, py, tmax) {
  const m = [1, 1, 1];
  const cx = px / 3, cy = py / 3; // centre of mass
  const x = [-0.5 - cx, -cy, 0, 0.5 - cx, -cy, 0, px - cx, py - cy, 0];
  const v = new Array(9).fill(0);
  const s = new NBody(m, x, v, { tol: TOL });
  let k = -1, tEsc = 0;
  try {
    while (s.t < tmax && s.steps < 150000) {
      s.step();
      if ((s.steps & 3) === 0) {
        const c = classify(s.x, s.v, m, 1, 3);
        if (c.state === "escape") { k = c.third; tEsc = s.t; break; }
      }
    }
    if (k === -1 && s.t < tmax) return { k: -2, t: s.t };
    let dE = s.relEnergyError();
    for (let extra = 0; extra < 64 && dE > 1e-6; extra++) { s.step(); dE = Math.min(dE, s.relEnergyError()); }
    if (dE > 1e-6) return { k: -2, t: s.t };
  } catch (e) {
    return { k: -2, t: s.t };
  }
  return { k, t: k >= 0 ? tEsc : s.t };
}

if (isMainThread) {
  const W = +(process.argv[2] || 400), H = +(process.argv[3] || 240), TMAX = +(process.argv[4] || 200);
  const nW = os.cpus().length;
  const half = new Uint8Array(W * H);
  const times = new Float32Array(W * H);
  let done = 0, next = 0;
  const t0 = Date.now();
  const rows = [...Array(H).keys()];
  const workers = [];
  const finish = () => {
    // row 0 = top (y = Y_MAX), column 0 = left (x = -X_MAX)
    const full = half;
    const counts = [0, 0, 0, 0, 0];
    for (const b of half) counts[b === 255 ? 3 : b === 254 ? 4 : b >> 6]++;
    const out = {
      width: W, height: H, xMin: -X_MAX, xMax: X_MAX, yMin: 0, yMax: Y_MAX,
      tMax: TMAX, tLo: T_LO, tol: TOL,
      note: "byte = escaper*64 + q (q = 63*ln(t/tLo)/ln(tMax/tLo)), 255 = no escape by tMax, 254 = unresolved (step cap, failure, or energy error > 1e-6)",
      counts: { body1: counts[0], body2: counts[1], body3: counts[2], none: counts[3], unresolved: counts[4] },
      data: Buffer.from(full).toString("base64"),
    };
    fs.writeFileSync(path.join(__dirname, "fatemap.json"), JSON.stringify(out));
    fs.writeFileSync(path.join(__dirname, "fatemap_times.bin"), Buffer.from(times.buffer));
    console.log(`done ${W}x${H} in ${((Date.now() - t0) / 1000).toFixed(1)}s`, out.counts);
  };
  for (let w = 0; w < nW; w++) {
    const worker = new Worker(__filename, { workerData: { W, H, TMAX } });
    const feed = () => {
      if (next < H) worker.postMessage(rows[next++]);
      else worker.postMessage(-1);
    };
    worker.on("message", (msg) => {
      for (let c = 0; c < W; c++) { half[msg.r * W + c] = msg.bytes[c]; times[msg.r * W + c] = msg.times[c]; }
      done++;
      if (done % Math.max(1, Math.floor(H / 10)) === 0)
        console.log(`${Math.round((100 * done) / H)}%  ${((Date.now() - t0) / 1000).toFixed(0)}s`);
      feed();
    });
    worker.on("exit", () => { workers.splice(workers.indexOf(worker), 1); if (!workers.length) finish(); });
    workers.push(worker);
    feed();
  }
} else {
  const { W, H, TMAX } = workerData;
  parentPort.on("message", (r) => {
    if (r < 0) return process.exit(0);
    const py = pixelY(r, H);
    const bytes = new Uint8Array(W), times = new Float32Array(W);
    for (let c = 0; c < W; c++) {
      const px = pixelX(c, W);
      const f = fate(px, py, TMAX);
      bytes[c] = f.k === -2 ? 254 : f.k < 0 ? 255 : f.k * 64 + quantize(f.t, TMAX);
      times[c] = f.k < 0 ? -1 : f.t;
    }
    parentPort.postMessage({ r, bytes, times });
  });
}
