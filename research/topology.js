/* Topology of a three-body orbit on Montgomery's shape sphere.
 *
 * Shape-sphere point (Li & Liao 2017, Šuvakov & Dmitrašinović 2013):
 *   ρ = (r1 - r2)/√2,  λ = (r1 + r2 - 2 r3)/√6,  R² = ρ² + λ²
 *   n = (2ρ·λ, λ² - ρ², 2(ρ×λ)_z) / R²
 * Binary collisions are three punctures on the equator (n_z = 0, the collinear configurations):
 *   C12 = (0, 1, 0), C13 = (-√3/2, -1/2, 0), C23 = (√3/2, -1/2, 0).
 * The equator splits into three arcs, labelled by the body in the middle of the line ("syzygy").
 *
 * Free-group word: project stereographically from one puncture P; the equator becomes a line
 * through the other two punctures, cut into L | M | R (M = the arc opposite P). Reading the
 * crossings of that line in (down, up) pairs gives letters:
 *   (R,M)=a  (L,M)=b  (M,R)=A  (M,L)=B  (L,R)=bA  (R,L)=aB  (X,X)=identity
 * with a/b the loops around the right/left punctures. Which puncture is P, which arc is L, and
 * which hemisphere is "up" are conventions; calibrate() picks the one matching published words.
 */
"use strict";

function shapePoint(x) {
  const r1 = [x[0], x[1]], r2 = [x[3], x[4]], r3 = [x[6], x[7]];
  const rho = [(r1[0] - r2[0]) / Math.SQRT2, (r1[1] - r2[1]) / Math.SQRT2];
  const s6 = Math.sqrt(6);
  const lam = [(r1[0] + r2[0] - 2 * r3[0]) / s6, (r1[1] + r2[1] - 2 * r3[1]) / s6];
  const R2 = rho[0] * rho[0] + rho[1] * rho[1] + lam[0] * lam[0] + lam[1] * lam[1];
  return [
    (2 * (rho[0] * lam[0] + rho[1] * lam[1])) / R2,
    (lam[0] * lam[0] + lam[1] * lam[1] - rho[0] * rho[0] - rho[1] * rho[1]) / R2,
    (2 * (rho[0] * lam[1] - rho[1] * lam[0])) / R2,
  ];
}

/* which body lies between the other two in a (nearly) collinear configuration */
function middleBody(x) {
  const d = (i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1]);
  const ex = [0, 1, 2].map((i) => {
    const [j, k] = [0, 1, 2].filter((q) => q !== i);
    return d(i, j) + d(i, k) - d(j, k);
  });
  return ex.indexOf(Math.min(...ex)) + 1; // 1, 2 or 3
}

/* Equator crossings over one period [t0, t0 + T]: [{t, body, dir}] with dir = -1 for
   n_z going + → − ("down") and +1 for − → +. Samples each solver step finely (Hermite) and
   refines each sign change by bisection on exact states (stateAt). */
function crossings(NBody, m, x0, v0, T, opts) {
  opts = opts || {};
  const tol = opts.tol ?? 1e-12, sub = opts.sub ?? 24;
  const sim = new NBody(m, x0, v0, { tol });
  // start away from t = 0: catalogue orbits begin exactly on the equator (a collinear configuration)
  const t0 = opts.t0 ?? 1e-3 * T;
  const out = [], buf = new Float64Array(9), st = new Float64Array(18);
  let prevT = null, prevZ = null;
  const nzAt = (t, exact) => {
    if (exact) { sim.stateAt(t, st); return shapePoint(st)[2]; }
    sim.sample(t, buf); return shapePoint(buf)[2];
  };
  while (sim.t < t0 + T) {
    sim.step();
    const a = sim.prev[18], b = Math.min(sim.t, t0 + T);
    const ts = [];
    for (let k = 1; k <= sub; k++) ts.push(a + ((sim.t - a) * k) / sub);
    if (sim.t >= t0 + T) ts.push(t0 + T); // make sure the window's end is sampled
    for (const t of ts.sort((p, q) => p - q)) {
      if (t < t0 || t > t0 + T) continue;
      const z = nzAt(t, false);
      if (prevZ !== null && Math.sign(z) !== Math.sign(prevZ) && z !== 0) {
        // bisect on exact states between prevT and t (both inside the current step when possible)
        let lo = Math.max(prevT, a), hi = t, zlo = nzAt(lo, true);
        for (let it = 0; it < 50 && hi - lo > 1e-12 * Math.max(1, hi); it++) {
          const mid = 0.5 * (lo + hi), zm = nzAt(mid, true);
          if (Math.sign(zm) === Math.sign(zlo)) { lo = mid; zlo = zm; } else hi = mid;
        }
        sim.stateAt(0.5 * (lo + hi), st);
        out.push({ t: 0.5 * (lo + hi), body: middleBody(st), dir: prevZ > 0 ? -1 : 1 });
      }
      prevT = t; prevZ = z;
    }
    if (b >= t0 + T) break;
  }
  return out;
}

/* ---- free-group words ---- */
const inv = (c) => (c === c.toLowerCase() ? c.toUpperCase() : c.toLowerCase());
function freeReduce(w) {
  const s = [];
  for (const c of w) { if (s.length && s[s.length - 1] === inv(c)) s.pop(); else s.push(c); }
  return s.join("");
}
function cyclicReduce(w) {
  w = freeReduce(w);
  while (w.length > 1 && w[0] === inv(w[w.length - 1])) w = w.slice(1, -1);
  return w;
}
const inverseWord = (w) => [...w].reverse().map(inv).join("");
/* canonical representative of the conjugacy class of w and of its inverse (time reversal) */
function canonical(w) {
  w = cyclicReduce(w);
  let best = null;
  for (const u of [w, inverseWord(w)])
    for (let i = 0; i < u.length; i++) {
      const r = u.slice(i) + u.slice(0, i);
      if (best === null || r < best) best = r;
    }
  return best ?? "";
}
/* expand catalogue notation such as "(ab)2ABBA(ba)2BAAB" or "b2(ABab)2A2" */
function expand(s) {
  s = s.replace(/\s+/g, "");
  let prev;
  do {
    prev = s;
    s = s.replace(/\(([abAB]+)\)(\d+)/g, (_, g, n) => g.repeat(+n));
    s = s.replace(/\(([abAB]+)\)/g, "$1");
    s = s.replace(/([abAB])(\d+)/g, (_, c, n) => c.repeat(+n));
  } while (s !== prev);
  return s;
}

/* word from crossings under a convention {pole: 1|2|3 (puncture C12, C13, C23), flip: bool
   (swap L/R), up: ±1 (which crossing direction counts as "down")} */
const PUNCTURES = { 1: [1, 2], 2: [1, 3], 3: [2, 3] };
function wordFrom(cr, conv) {
  const [p, q] = PUNCTURES[conv.pole];
  const mid = [1, 2, 3].find((b) => b !== p && b !== q); // arc opposite the pole: that body in the middle
  // the two arcs touching the pole: body p in the middle, or body q in the middle
  const seg = (b) => (b === mid ? "M" : (b === p) !== conv.flip ? "L" : "R");
  // rotate so the sequence starts with a "down" crossing
  const down = -conv.up;
  let k = cr.findIndex((c) => c.dir === down);
  if (k < 0 || cr.length % 2) return null;
  const seq = cr.slice(k).concat(cr.slice(0, k));
  const LET = { RM: "a", LM: "b", MR: "A", ML: "B", LR: "bA", RL: "aB" };
  let w = "";
  for (let i = 0; i < seq.length; i += 2) {
    const key = seg(seq[i].body) + seg(seq[i + 1].body);
    if (key[0] !== key[1]) w += LET[key];
  }
  return cyclicReduce(w);
}
const CONVENTIONS = [];
for (const pole of [1, 2, 3]) for (const flip of [false, true]) for (const up of [1, -1]) CONVENTIONS.push({ pole, flip, up });

module.exports = { shapePoint, middleBody, crossings, freeReduce, cyclicReduce, inverseWord, canonical, expand, wordFrom, CONVENTIONS };
