/* Data for the findings page (research/web/dist/ladder.html).
 * Integrates every orbit shown on the page over one period, records a display path (decimated,
 * quantized), the syzygy letters with their times, closure and stability; gathers continuation
 * branches, the search results per mass, the published catalogue (for the discovery map) and the
 * return-distance maps. Usage: node build_data.js [out.json] */
"use strict";
const fs = require("fs"), path = require("path");
const R = path.join(__dirname, "..");
const { NBody } = require(path.join(R, "../web/solver.js"));
const S = require(path.join(R, "slice.js")), TP = require(path.join(R, "topology.js")), RF = require(path.join(R, "refine.js"));
const CONV = { pole: 1, flip: false, up: 1 };
const out = process.argv[2] || path.join(__dirname, "ladder_data.json");
const readJSON = (f) => JSON.parse(fs.readFileSync(path.join(R, f), "utf8"));
const readJSONL = (f) => fs.existsSync(path.join(R, f)) ? fs.readFileSync(path.join(R, f), "utf8").trim().split("\n").filter(Boolean).map((l) => JSON.parse(l)) : [];
const b64 = (typed) => Buffer.from(typed.buffer, typed.byteOffset, typed.byteLength).toString("base64");

/* ---------- orbits shown on the page ---------- */
const claims = readJSON("results/claims.json");
const stab = Object.fromEntries(readJSON("results/claims_stability.json").map((s) => [s.id, s]));
const orbits = [];
for (const c of claims) {
  const known = c.status.startsWith("control");
  let id = c.id, k = c.k;
  if (id.startsWith("catalogue")) { k = c.word.length / 4 - 1; id = `L${k}@0.5`; }
  orbits.push({ id, m3: c.m3, k, v1: c.v1, v2: c.v2, T: c.T, status: known ? "known" : "new", note: c.status, stabId: c.id,
    match: known ? (c.id.startsWith("catalogue") ? c.id.replace("catalogue ", "Li, Jing & Liao 2018, ") : "word catalogued (I.A family)") : null, extra: id.startsWith("X") });
}
// more rungs found later by the grid at m3 = 0.15, 0.2, 0.3 or by continuation (results/claims_extra.json,
// closure checked by the same independent integrator), and catalogued-word rungs at those masses
const stabExtra = fs.existsSync(path.join(R, "results/claims_extra_stability.json")) ? readJSON("results/claims_extra_stability.json") : [];
for (const s of stabExtra) stab[s.id] = s;
if (fs.existsSync(path.join(R, "results/claims_extra.json")))
  for (const c of readJSON("results/claims_extra.json")) orbits.push({ id: c.id, m3: c.m3, k: c.k, v1: c.v1, v2: c.v2, T: c.T, status: "new", note: c.status, stabId: c.id, match: null, later: true });
if (fs.existsSync(path.join(R, "results/rungs_known_extra.json")))
  for (const c of readJSON("results/rungs_known_extra.json")) orbits.push({ id: c.id, m3: c.m3, k: c.k, v1: c.v1, v2: c.v2, T: c.T, status: "known", note: c.status, match: "word catalogued (I.A family)", firstAtMass: true });
// the k = 2 rung at m3 = 0.1 comes from the m3 = 0.1 sweep (its word is catalogued at m3 = 0.5)
for (const o of readJSONL("runs/m01_N800_T100_orbits.jsonl"))
  if (o.canon === "ABaBaBabAbAb") orbits.push({ id: "L2@0.1", m3: 0.1, k: 2, v1: o.v1, v2: o.v2, T: o.T, status: "known", note: "sweep at m3 = 0.1", match: "word catalogued at m3 = 0.5 (I.A i.c. 3)" });
// rungs 1 and 2 at m3 = 0.75 and rung 1 at m3 = 1, from continuation (results/anchors.json)
for (const a of readJSON("results/anchors.json"))
  orbits.push({ id: a.id, m3: a.m3, k: a.k, v1: a.v1, v2: a.v2, T: a.T, status: "known", note: a.via,
    match: a.published ? a.published.source : a.note, anchor: a.published ? { v1: a.published.v1, v2: a.published.v2, T: a.published.T } : null, firstAtMass: !a.published });

function letters(cr) { // syzygy letters with the time each one is completed (same rule as TP.wordFrom)
  const [p, q] = [1, 2], mid = 3, seg = (b) => (b === mid ? "M" : (b === p) !== CONV.flip ? "L" : "R");
  const k = cr.findIndex((c) => c.dir === -CONV.up);
  const seq = cr.slice(k).concat(cr.slice(0, k));
  const LET = { RM: "a", LM: "b", MR: "A", ML: "B", LR: "bA", RL: "aB" }, outL = [];
  for (let i = 0; i < seq.length; i += 2) {
    const key = seg(seq[i].body) + seg(seq[i + 1].body);
    outL.push([+seq[i + 1].t.toFixed(5), key[0] !== key[1] ? LET[key] : ""]);
  }
  return { letters: outL.sort((a, b) => a[0] - b[0]), syzygies: cr.map((c) => [+c.t.toFixed(5), c.body]) };
}

function track(m3, v1, v2, T) {
  const { m, x, v } = S.ic(m3, v1, v2);
  const sim = new NBody(m, x, v, { tol: 1e-12 });
  const ts = [0], ps = [[...x.filter((_, i) => i % 3 !== 2)]];
  const buf = new Float64Array(9), sub = 16;
  while (sim.t < T) {
    sim.step();
    const a = sim.prev[18], b = Math.min(sim.t, T);
    for (let j = 1; j <= sub; j++) {
      const t = a + ((b - a) * j) / sub;
      sim.sample(t, buf);
      ts.push(t); ps.push([buf[0], buf[1], buf[3], buf[4], buf[6], buf[7]]);
    }
  }
  // closure, from the exact state at T
  const st = new Float64Array(18); sim.stateAt(T, st);
  const y0 = [...x, ...v]; let cl = 0;
  for (let i = 0; i < 18; i++) cl += (st[i] - y0[i]) ** 2;
  // greedy decimation: keep a sample whenever linear interpolation in time would err by > eps
  const eps = 2.5e-4, keep = [0];
  let i = 0;
  while (i < ts.length - 1) {
    let j = i + 2;
    for (; j < ts.length; j++) {
      let bad = false;
      for (let q = i + 1; q < j && !bad; q++) {
        const f = (ts[q] - ts[i]) / (ts[j] - ts[i]);
        for (let c = 0; c < 6; c += 2) {
          const ex = ps[i][c] + f * (ps[j][c] - ps[i][c]) - ps[q][c], ey = ps[i][c + 1] + f * (ps[j][c + 1] - ps[i][c + 1]) - ps[q][c + 1];
          if (ex * ex + ey * ey > eps * eps) { bad = true; break; }
        }
      }
      if (bad || j - i > 40) break;
    }
    i = Math.max(i + 1, j - 1); keep.push(i);
  }
  let lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
  for (const p of ps) for (let c = 0; c < 6; c++) { lo[c & 1] = Math.min(lo[c & 1], p[c]); hi[c & 1] = Math.max(hi[c & 1], p[c]); }
  const span = Math.max(hi[0] - lo[0], hi[1] - lo[1]), cx = (lo[0] + hi[0]) / 2, cy = (lo[1] + hi[1]) / 2;
  const tq = new Float32Array(keep.length), xq = new Int16Array(keep.length * 6);
  keep.forEach((s, n) => {
    tq[n] = ts[s];
    for (let c = 0; c < 6; c++) xq[n * 6 + c] = Math.round(((ps[s][c] - (c & 1 ? cy : cx)) / span) * 2 * 32760);
  });
  return { n: keep.length, closure: Math.sqrt(cl), box: [cx, cy, span], t: b64(tq), xy: b64(xq), raw: ts.length };
}

for (const o of orbits) if (o.status === "known") o.firstAtMass = !(o.match || "").startsWith("Li");
for (const o of orbits) {
  // tighten the published 10-digit values (and anything not refined at 1e-13) before tracing
  if (o.status === "known") {
    const f = RF.refine(o.m3, [o.v1, o.v2, o.T], { maxIt: 20 });
    if (f.ok) [o.v1, o.v2, o.T] = f.p;
  }
  o.E = S.energy(o.m3, o.v1, o.v2); o.Tstar = o.T * Math.pow(-o.E, 1.5);
  const { m, x, v } = S.ic(o.m3, o.v1, o.v2);
  const cr = TP.crossings(NBody, m, x, v, o.T, { tol: 1e-13 });
  o.word = TP.wordFrom(cr, CONV); o.canon = TP.canonical(o.word);
  Object.assign(o, letters(cr));
  const raw = o.letters.map((l) => l[1]).join("");
  o.wordRawMatches = TP.canonical(raw) === o.canon && raw.length === o.canon.length;
  // the word as the page spells it: one entry per letter, with the time it is written, in time
  // order from t = 0 (the starting line-up). For the ladder this reads (Ba)^(k+1) (bA)^(k+1).
  const flat = []; for (const [t, s] of o.letters) for (const ch of s) flat.push([t, ch]);
  o.spell = flat;
  o.spellCanonical = TP.canonical(flat.map((f) => f[1]).join("")) === o.canon;
  delete o.letters;
  Object.assign(o, track(o.m3, o.v1, o.v2, o.T));
  const s = stab[o.stabId];
  if (s) { o.maxMult = s.max_multiplier; o.lyap = s.lyapunov_per_time; o.mults = s.multipliers.map((z) => [+z[0].toPrecision(6), +z[1].toPrecision(6)]); }
  delete o.stabId;
  console.log(`${o.id.padEnd(9)} k=${o.k} T=${o.T.toFixed(4)} word ${o.canon} raw-ok ${o.wordRawMatches} spell-canon ${o.spellCanonical} samples ${o.n}/${o.raw} closure ${o.closure.toExponential(2)} maxMult ${o.maxMult ? o.maxMult.toFixed(2) : "-"}`);
}

/* ---------- continuation branches (m3, T, T*) ---------- */
const branches = [];
for (let k = 1; k <= 8; k++) {
  const files = [`runs/contT_k${k}.json`, ...(k === 3 ? ["runs/contT_k3_hi_rev.json", "runs/contT_k3_hi.json"] : [])].filter((f) => fs.existsSync(path.join(R, f)));
  let pts = [];
  for (const f of files) pts.push(...readJSON(f));
  for (const extra of [`runs/contT_k${k}_down.json`, `runs/contM_k${k}.json`]) if (fs.existsSync(path.join(R, extra))) pts.push(...readJSON(extra));
  pts = pts.filter((p) => p.m3 > 0 && S.energy(p.m3, p.v1, p.v2) < 0).sort((a, b) => a.T - b.T);
  const seen = new Set();
  pts = pts.filter((p) => { const key = p.T.toFixed(3); if (seen.has(key)) return false; seen.add(key); return true; });
  if (pts.length < 2) continue;
  const step = Math.max(1, Math.floor(pts.length / 160));
  const pick = pts.filter((_, i) => i % step === 0 || i === pts.length - 1);
  branches.push({ k, pts: pick.map((p) => { const E = S.energy(p.m3, p.v1, p.v2); return [+p.m3.toFixed(5), +p.T.toFixed(4), +(p.T * Math.pow(-E, 1.5)).toFixed(4)]; }) });
}

/* ---------- published catalogue (discovery map background) ---------- */
const catN = readJSON("data/catalog_numeric.json");
const catRows = catN.map((r) => [r.m3, +(+r.Tstar).toFixed(3), r.Lf, r.label]).filter((r) => isFinite(r[1]));

/* ---------- sweeps ---------- */
const sweeps = [];
for (const f of fs.readdirSync(path.join(R, "runs")).filter((f) => /^m\d+_N800_T100_orbits\.jsonl$/.test(f)).sort()) {
  const rows = readJSONL("runs/" + f);
  if (!rows.length) continue;
  const m3 = rows[0].m3;
  const log = fs.existsSync(path.join(R, "runs/" + f.replace("_orbits.jsonl", ".log"))) ? fs.readFileSync(path.join(R, "runs/" + f.replace("_orbits.jsonl", ".log")), "utf8") : "";
  const last = (re) => { const all = [...log.matchAll(re)]; return all.length ? all[all.length - 1][1] : null; };
  const cand = last(/(\d+) candidate rows/g), jobs = last(/(\d+) refinement jobs/g); // logs can hold reruns: use the last
  const converged = last(/done: \d+ jobs, (\d+) converged/g);
  // drop repeated traversals (u^k found at k T): keep primitive orbits only
  const prim = rows.filter((o) => {
    if (o.power === 1) return true;
    const r = RF.residual(m3, [o.v1, o.v2, o.T / o.power], 1e-12);
    return !(r && Math.hypot(...r) < 1e-7);
  });
  const ladderWords = new Set(orbits.map((o) => o.canon));
  sweeps.push({ m3, candidates: cand ? +cand : null, jobs: jobs ? +jobs : null, converged: converged ? +converged : null,
    orbits: prim.map((o) => ({ v1: +o.v1.toFixed(6), v2: +o.v2.toFixed(6), T: +o.T.toFixed(4), Ts: +o.Tstar.toFixed(3), Lf: o.Lf, known: o.known, ladder: ladderWords.has(o.canon), canon: o.canon, match: o.matches[0] || null, minSep: +o.minSep.toPrecision(3), power: o.power })) });
  console.log(`sweep m3=${m3}: ${rows.length} orbits, ${prim.length} primitive, ${prim.filter((o) => !o.known).length} not in catalogue`);
}

/* ---------- return-distance maps ---------- */
const maps = [];
for (const f of ["runs/map_m025_N400.json", "runs/map_m05_N400.json"]) {
  if (!fs.existsSync(path.join(R, f))) { console.warn(`missing ${f} (run map.js); the page will have no map for it`); continue; }
  const M = readJSON(f), N = M.N;
  const d = new Uint8Array(N * N), T = new Uint8Array(N * N);
  for (let i = 0; i < N * N; i++) {
    if (M.d[i] < 0) { d[i] = 0; continue; } // unbound
    const L = Math.log10(Math.max(M.d[i], 1e-3)); // 1e-3 .. 10
    d[i] = 1 + Math.max(0, Math.min(253, Math.round(((L + 3) / 4) * 253)));
    T[i] = M.fate[i] === 1 ? 255 : Math.round((M.T[i] / M.Tmax) * 254);
  }
  maps.push({ m3: M.m3, N, Tmax: M.Tmax, vmax: M.vmax, d: b64(d), T: b64(T) });
}

/* ---------- verification summary and scan size ---------- */
const verify = fs.existsSync(path.join(R, "results/verification.json")) ? readJSON("results/verification.json").cards : null;
const bound = (m3, N) => { const vm = S.vmax(m3); let n = 0; for (let r = 0; r < N; r++) for (let c = 0; c < N; c++) if (S.energy(m3, (vm * (c + 0.5)) / N, (vm * (r + 0.5)) / N) < 0) n++; return n; };
const scan = { perMass: sweeps.map((s) => ({ m3: s.m3, points: bound(s.m3, 800) })) };
scan.total = scan.perMass.reduce((a, b) => a + b.points, 0);

const data = { built: new Date().toISOString().slice(0, 10), orbits, branches, catalogue: catRows, sweeps, maps, verify, scan };
fs.writeFileSync(out, JSON.stringify(data));
console.log(`wrote ${out} (${(fs.statSync(out).size / 1024).toFixed(0)} KB)`);
