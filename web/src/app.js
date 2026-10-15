"use strict";
(function () {
  /* ================================================================ constants */
  const BODY = ["#c2841a", "#3a91d4", "#de5577"]; // validated chart steps
  const GLOW = ["#f6c66a", "#8fd0ff", "#ff9bb3"]; // star cores
  const GLOW_RGB = [[246, 198, 106], [143, 208, 255], [255, 155, 179]];
  const $ = (id) => document.getElementById(id);
  const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
  const reduceMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const DAY_PER_TU = 365.25 / (2 * Math.PI); // 1 AU, 1 M_sun, G = 1  →  t* = yr / 2π

  function fmt(x, d = 3) {
    if (!isFinite(x)) return x > 0 ? "∞" : "–";
    const ax = Math.abs(x);
    if (ax !== 0 && (ax < 1e-3 || ax >= 1e5)) return x.toExponential(Math.max(0, d - 1)).replace("e", "e");
    return x.toFixed(d);
  }
  function sci(x) {
    if (!isFinite(x) || x <= 0) return x === 0 ? "0" : "–";
    const e = Math.floor(Math.log10(x));
    const m = x / Math.pow(10, e);
    return `${m.toFixed(1)}e${e}`;
  }
  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6d2b79f5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* ================================================================ initial conditions */
  function toCom(m, x, v) {
    const M = m[0] + m[1] + m[2];
    for (let k = 0; k < 3; k++) {
      let X = 0, V = 0;
      for (let i = 0; i < 3; i++) { X += m[i] * x[3 * i + k]; V += m[i] * v[3 * i + k]; }
      for (let i = 0; i < 3; i++) { x[3 * i + k] -= X / M; v[3 * i + k] -= V / M; }
    }
    return { m: m.slice(), x, v };
  }
  function keplerPair(m1, m2, a, e) { // pericentre, xy plane
    const M = m1 + m2, r = a * (1 - e), vr = Math.sqrt((M * (1 + e)) / (a * (1 - e)));
    return [[(-m2 / M) * r, 0, 0], [0, (-m2 / M) * vr, 0], [(m1 / M) * r, 0, 0], [0, (m1 / M) * vr, 0]];
  }
  function lagrangeIC(m, a = 1, e = 0) {
    const x = [], v = [];
    const M = m[0] + m[1] + m[2];
    const ang = [90, 210, 330].map((d) => (d * Math.PI) / 180);
    const P = ang.map((t) => [(a / Math.sqrt(3)) * Math.cos(t), (a / Math.sqrt(3)) * Math.sin(t)]);
    const cx = (m[0] * P[0][0] + m[1] * P[1][0] + m[2] * P[2][0]) / M;
    const cy = (m[0] * P[0][1] + m[1] * P[1][1] + m[2] * P[2][1]) / M;
    const w = Math.sqrt(M / (a * a * a)) * Math.sqrt(1 + e);
    const A = a / (1 - e); // semi-major axis of the homographic Kepler motion
    for (let i = 0; i < 3; i++) {
      const px = P[i][0] - cx, py = P[i][1] - cy;
      x.push(px, py, 0); v.push(-w * py, w * px, 0);
    }
    const r = toCom(m, x, v);
    r.period = 2 * Math.PI * Math.sqrt((A * A * A) / M);
    return r;
  }
  function eulerIC(m, a = 2) {
    const [m1, m2, m3] = m;
    const acc = (l) => [m2 + m3 / (1 + l) ** 2, -m1 + m3 / (l * l), -m1 / (1 + l) ** 2 - m2 / (l * l)];
    const f = (l) => { const [a1, a2, a3] = acc(l); return a2 - a1 - (a3 - a2) / l; };
    let lo = 1e-6, hi = 1e6, flo = f(lo);
    for (let k = 0; k < 200; k++) {
      const mid = Math.sqrt(lo * hi), fm = f(mid);
      if (Math.sign(fm) === Math.sign(flo)) { lo = mid; flo = fm; } else hi = mid;
    }
    const lam = Math.sqrt(lo * hi);
    const [a1, a2] = acc(lam);
    const s = a / (1 + lam);
    const om = Math.sqrt(-(a2 - a1) / (s * s * s));
    const xs = [0, s, s * (1 + lam)];
    const xc = (m1 * xs[0] + m2 * xs[1] + m3 * xs[2]) / (m1 + m2 + m3);
    const x = [], v = [];
    for (let i = 0; i < 3; i++) { x.push(xs[i] - xc, 0, 0); v.push(0, om * (xs[i] - xc), 0); }
    const r = toCom(m, x, v);
    r.lambda = lam;
    r.period = (2 * Math.PI) / om;
    return r;
  }
  function suvakovIC(p1, p2) {
    return toCom([1, 1, 1], [-1, 0, 0, 1, 0, 0, 0, 0, 0], [p1, p2, 0, p1, p2, 0, -2 * p1, -2 * p2, 0]);
  }
  function freeFallIC(px, py) { // identical arithmetic to build_fatemap.js
    const cx = px / 3, cy = py / 3;
    return { m: [1, 1, 1], x: [-0.5 - cx, -cy, 0, 0.5 - cx, -cy, 0, px - cx, py - cy, 0], v: [0, 0, 0, 0, 0, 0, 0, 0, 0] };
  }
  function randomIC(seed) {
    const r = mulberry32(seed);
    const m = [0.5 + r(), 0.5 + r(), 0.5 + r()];
    const x = [], v = [];
    for (let i = 0; i < 3; i++) {
      const ang = 2 * Math.PI * r(), rad = 0.4 + 0.8 * Math.sqrt(r());
      x.push(rad * Math.cos(ang), rad * Math.sin(ang), 0);
      const va = 2 * Math.PI * r(), vm = r();
      v.push(vm * Math.cos(va), vm * Math.sin(va), 0);
    }
    const s = toCom(m, x, v);
    // rescale velocities to a virial ratio Q = T/U = 0.25 (bound, chaotic)
    let T = 0, U = 0;
    for (let i = 0; i < 3; i++) T += 0.5 * m[i] * (s.v[3 * i] ** 2 + s.v[3 * i + 1] ** 2);
    for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++)
      U += (m[i] * m[j]) / Math.hypot(s.x[3 * i] - s.x[3 * j], s.x[3 * i + 1] - s.x[3 * j + 1]);
    const f = Math.sqrt((0.25 * U) / Math.max(T, 1e-12));
    for (let k = 0; k < 9; k++) s.v[k] *= f;
    return s;
  }
  function routh(m) {
    const lhs = 27 * (m[0] * m[1] + m[1] * m[2] + m[2] * m[0]), rhs = (m[0] + m[1] + m[2]) ** 2;
    return { lhs, rhs, stable: lhs < rhs };
  }

  const orb = (k) => (ORBITS && ORBITS[k]) || null;
  const PRESETS = [];
  function addSuvakov(id, name, key, blurb, speed) {
    const o = orb(key);
    if (!o) return;
    PRESETS.push({
      id, name, group: "Periodic orbits", masses: [1, 1, 1], massMode: "keep",
      make: () => suvakovIC(o.p1, o.p2), period: o.T, speed, trail: Math.min(o.T, 30),
      blurb: `${blurb} Period T = ${o.T.toFixed(6)}. Starting values refined by Newton shooting so the orbit closes to ${sci(o.residual)}.`,
      source: "Šuvakov & Dmitrašinović 2013, PRL 110, 114301 (refined)",
    });
  }
  PRESETS.push({
    id: "figure8", name: "Figure-eight", group: "Periodic orbits", masses: [1, 1, 1], massMode: "keep",
    make: () => { const x1 = [-0.97000436, 0.24308753], v3 = [-0.93240737, -0.86473146];
      return toCom([1, 1, 1], [x1[0], x1[1], 0, -x1[0], -x1[1], 0, 0, 0, 0], [-v3[0] / 2, -v3[1] / 2, 0, -v3[0] / 2, -v3[1] / 2, 0, v3[0], v3[1], 0]); },
    period: 6.32591398, speed: 1.4, trail: 6.33, preroll: true,
    blurb: "Three equal masses chase one another around a single figure-eight. It is <b>linearly stable</b>: nudge it and it keeps its shape. Try the butterfly twin and compare with the Pythagorean problem.",
    source: "Moore 1993; Chenciner & Montgomery 2000",
  });
  PRESETS.push({
    id: "lagrange", name: "Lagrange triangle", group: "Exact solutions (any masses)", masses: [1, 1, 1], massMode: "regen",
    make: (m) => lagrangeIC(m, 1, 0), speed: 1.2, trail: 9,
    blurb: "An equilateral triangle that rotates rigidly. It is an exact solution <b>for any masses</b>: move the sliders and the page rebuilds it. ",
    routh: true, source: "Lagrange 1772; stability: Routh 1875",
  });
  PRESETS.push({
    id: "lagrangeE", name: "Lagrange, pulsating", group: "Exact solutions (any masses)", masses: [1, 0.5, 0.25], massMode: "regen",
    make: (m) => lagrangeIC(m, 1, 0.45), speed: 1.2, trail: 12,
    blurb: "The same triangle, but each body follows a Kepler ellipse (e = 0.45) about the barycentre. The triangle grows and shrinks while staying equilateral. ",
    routh: true, source: "Lagrange 1772: homographic solutions",
  });
  PRESETS.push({
    id: "euler", name: "Euler line", group: "Exact solutions (any masses)", masses: [1, 2, 1.5], massMode: "regen",
    make: (m) => eulerIC(m, 2), speed: 1.2, trail: 8,
    blurb: "The bodies sit on a straight line that rotates rigidly. The spacing ratio λ is the positive root of Euler's quintic, so it changes as you move the mass sliders. Collinear solutions are always unstable, so it eventually falls apart.",
    euler: true, source: "Euler 1767",
  });
  addSuvakov("butterfly1", "Butterfly I", "butterfly_I", "Two wings traced back and forth.", 1.4);
  addSuvakov("butterfly2", "Butterfly II", "butterfly_II", "A wider-winged sibling of Butterfly I.", 1.4);
  addSuvakov("butterfly3", "Butterfly III", "butterfly_III", "Butterfly with nested loops.", 2);
  addSuvakov("moth1", "Moth I", "moth_I", "Interlocking loops with three-fold winding.", 2);
  addSuvakov("moth2", "Moth II", "moth_II", "A longer moth with doubled loops.", 3);
  addSuvakov("moth3", "Moth III", "moth_III", "Moth with a richer braid.", 3);
  addSuvakov("yinyang1a", "Yin-yang Ia", "yin_yang_Ia", "Two lobes that swap roles each half-period.", 2.5);
  addSuvakov("yinyang1b", "Yin-yang Ib", "yin_yang_Ib", "A tighter yin-yang.", 2);
  addSuvakov("goggles", "Goggles", "goggles", "Two lenses joined by a bridge.", 1.6);
  addSuvakov("dragonfly", "Dragonfly", "dragonfly", "Long thin wings swept by fast passes.", 3);
  addSuvakov("yarn", "Yarn", "yarn", "A long, tangled ball of three threads.", 5);
  addSuvakov("bumblebee", "Bumblebee", "bumblebee", "One of the longest-period orbits in the family.", 6);
  PRESETS.push({
    id: "pythagorean", name: "Pythagorean problem", group: "Astrophysics and chaos", masses: [3, 4, 5], massMode: "keep",
    make: () => toCom([3, 4, 5], [1, 3, 0, -2, -1, 0, 1, -1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]),
    speed: 2, trail: 7,
    blurb: "Masses 3, 4 and 5 start at rest on a 3-4-5 right triangle. After dozens of close encounters (one within 0.0006) the lightest body is ejected around t ≈ 60 and the other two leave as a tight binary.",
    source: "Burrau 1913; solved by Szebehely & Peters 1967",
  });
  PRESETS.push({
    id: "kepler16", name: "Kepler-16: two suns", group: "Astrophysics and chaos", masses: [0.6897, 0.20255, 3.18e-4], massMode: "keep",
    make: () => {
      const [xA, vA, xB, vB] = keplerPair(0.6897, 0.20255, 0.22431, 0.15944);
      const ap = 0.7048, vp = Math.sqrt((0.6897 + 0.20255 + 3.18e-4) / ap);
      return toCom([0.6897, 0.20255, 3.18e-4], [...xA, ...xB, 0, ap, 0], [...vA, ...vB, -vp, 0, 0]);
    },
    speed: 1.2, trail: 4.2, units: "days", scaleNote: "AU, M☉",
    blurb: "A real system: a Saturn-mass planet orbiting a pair of stars every 229 days, 245 light-years away. The binary's 41-day dance makes the planet's orbit wobble. Units: AU and solar masses; the clock reads days.",
    source: "Doyle et al. 2011, Science 333 (masses, orbits)",
  });
  PRESETS.push({
    id: "moon", name: "Star, planet and moon", group: "Astrophysics and chaos", masses: [1, 0.003, 3e-5], massMode: "keep",
    make: () => {
      const vp = Math.sqrt(1.003), am = 0.04, vm = Math.sqrt(0.00303 / am);
      return toCom([1, 0.003, 3e-5], [0, 0, 0, 1, 0, 0, 1 + am, 0, 0], [0, 0, 0, 0, vp, 0, 0, vp + vm, 0]);
    },
    speed: 1, trail: 6.3,
    blurb: "A hierarchical triple: the moon circles the planet well inside its Hill sphere, so the star barely disturbs it. Seen from the star, the moon's path is a looping epicycle.",
    source: "Hill stability: r_H = a (m / 3M)^1/3",
  });
  PRESETS.push({
    id: "exchange", name: "Binary–single exchange", group: "Astrophysics and chaos", masses: [1, 1, 1.2], massMode: "keep",
    make: () => {
      const phi = Math.PI / 12, b = 0.6, v0 = 0.4, D = 12, vc = Math.SQRT2 / 2, c = Math.cos(phi), s = Math.sin(phi);
      return toCom([1, 1, 1.2], [-0.5 * c, -0.5 * s, 0, 0.5 * c, 0.5 * s, 0, -D, b, 0], [vc * s, -vc * c, 0, -vc * s, vc * c, 0, v0, 0, 0]);
    },
    speed: 3.5, trail: 14,
    blurb: "A heavier star falls onto a circular binary. After a resonant scramble it captures a partner and throws one of the original pair out. Exchanges like this build compact binaries in dense star clusters.",
    source: "Hills 1975; Heggie 1975",
  });
  PRESETS.push({
    id: "kozai", name: "Kozai–Lidov cycles", group: "Astrophysics and chaos", masses: [1, 0.001, 1], massMode: "keep", tilt: 62,
    make: () => {
      const [x1, v1, x2, v2] = keplerPair(1, 0.001, 1, 0.05);
      const ao = 5, vo = Math.sqrt(2.001 / ao), i = (80 * Math.PI) / 180;
      return toCom([1, 0.001, 1], [...x1, ...x2, ao, 0, 0], [...v1, ...v2, 0, vo * Math.cos(i), vo * Math.sin(i)]);
    },
    speed: 28, trail: 40,
    blurb: "A planet around a star, with a companion star on an orbit tilted 80°. The planet's orbit trades inclination for eccentricity: watch it stretch into a needle and back. Fully 3D; drag to orbit the view.",
    source: "von Zeipel 1910; Lidov 1962; Kozai 1962",
  });
  PRESETS.push({
    id: "sitnikov", name: "Sitnikov problem", group: "Astrophysics and chaos", masses: [0.5, 0.5, 1e-5], massMode: "keep", tilt: 66,
    make: () => {
      const [x1, v1, x2, v2] = keplerPair(0.5, 0.5, 1, 0.2);
      return toCom([0.5, 0.5, 1e-5], [...x1, ...x2, 0, 0, 1.1], [...v1, ...v2, 0, 0, 0]);
    },
    speed: 2.5, trail: 30,
    blurb: "Two equal stars on an eccentric orbit; a light third body bobs up and down through their centre. The bobbing is fully chaotic: one of the simplest systems proven to be chaotic.",
    source: "Sitnikov 1960; Moser 1973",
  });
  PRESETS.push({
    id: "trojan", name: "Trojan at L4", group: "Astrophysics and chaos", masses: [0.999, 0.001, 1e-9], massMode: "keep", rotating: true,
    make: () => {
      const mu = 0.001, a = (78 * Math.PI) / 180, L4 = [-mu + Math.cos(a), Math.sin(a), 0]; // 18° ahead of L4
      const vv = (p) => [-p[1], p[0], 0];
      return toCom([1 - mu, mu, 1e-9], [-mu, 0, 0, 1 - mu, 0, 0, ...L4], [...vv([-mu, 0]), ...vv([1 - mu, 0]), ...vv(L4)]);
    },
    speed: 9, trail: 160,
    blurb: "A Sun–Jupiter pair and an asteroid released 78° ahead of the planet, 18° past the L4 point. In the co-rotating frame the gravity field becomes the effective potential, and the asteroid swings back and forth between about 46° and 78° in a tadpole around L4, every 76 time units. Lagrange's triangle again, made stable by Routh's criterion.",
    source: "Lagrange 1772; Jupiter's Trojans, Wolf 1906",
  });
  PRESETS.push({
    id: "random", name: "Random bound triple", group: "Astrophysics and chaos", masses: null, massMode: "keep",
    make: () => randomIC(state.seed), speed: 1.5, trail: 7, random: true,
    blurb: "Random masses and positions, velocities scaled so the system is bound (virial ratio 0.25). Almost every such system ends with an ejection.",
    source: "",
  });
  PRESETS.push({
    id: "freefall", name: "Free fall (from fate map)", group: "Astrophysics and chaos", masses: [1, 1, 1], massMode: "keep",
    make: () => freeFallIC(state.ff[0], state.ff[1]), speed: 2, trail: 8, tol: 1e-11, rMin: 3,
    blurb: "Equal masses released from rest. Pick a starting point for body 3 on the fate map below the simulator.",
    source: "Agekyan & Anosova 1967",
  });

  /* ================================================================ state */
  const state = {
    preset: null, base: null, masses: [1, 1, 1], custom: false, seed: 7, ff: [0.2, 0.6],
    sim: null, twin: null, playing: !reduceMotion, speed: 1, tDisp: 0, lag: false,
    opts: { trails: true, field: true, twin: false, rotating: false, velocity: false, forces: false, labels: true, autofit: true },
    tilt: 0, yaw: 0, view3d: false, trailDur: 6, tol: 1e-11, eps: 0, delta: 1e-9,
    trails: [], twinTrails: [], diag: null, sphere: [], periods: [], status: null, classifyEvery: 0,
    disp: new Float64Array(9), dispV: new Float64Array(9), twinDisp: new Float64Array(9), theta: 0, omega: 0,
    y0: null, sizeScale: 1, editing: null,
  };

  /* ================================================================ trails */
  class Trail {
    constructor(cap) { this.cap = cap; this.x = new Float32Array(cap); this.y = new Float32Array(cap); this.z = new Float32Array(cap);
      this.c = new Float32Array(cap); this.s = new Float32Array(cap); this.t = new Float64Array(cap); this.head = 0; this.len = 0; }
    push(x, y, z, t, c, s) { const i = this.head; this.x[i] = x; this.y[i] = y; this.z[i] = z; this.t[i] = t; this.c[i] = c; this.s[i] = s;
      this.head = (i + 1) % this.cap; if (this.len < this.cap) this.len++; }
    clear() { this.head = 0; this.len = 0; }
  }

  /* ================================================================ diagnostics */
  class Series {
    constructor(n) { this.n = n; this.t = []; this.y = []; for (let i = 0; i < n; i++) this.y.push([]); }
    push(t, vals) { this.t.push(t); for (let i = 0; i < this.n; i++) this.y[i].push(vals[i]);
      if (this.t.length > 2400) { this.t = this.t.filter((_, k) => k % 2 === 0); this.y = this.y.map((a) => a.filter((_, k) => k % 2 === 0)); } }
  }
  function newDiag() {
    return { E: new Series(1), R: new Series(3), S: new Series(1), lastT: -Infinity, lyap: { n: 0, st: 0, sy: 0, stt: 0, sty: 0, horizon: null }, minR: Infinity };
  }

  /* ================================================================ simulation */
  function currentPreset() { return state.preset; }

  function loadPreset(id, keepMasses) {
    const p = PRESETS.find((q) => q.id === id) || PRESETS[0];
    state.preset = p;
    state.custom = false;
    if (!keepMasses) state.masses = p.masses ? p.masses.slice() : null;
    let ic = p.massMode === "regen" ? p.make(state.masses) : p.make();
    if (!state.masses) state.masses = ic.m.slice();
    if (p.massMode === "keep" && keepMasses && state.masses) ic = toCom(state.masses.slice(), ic.x.slice(), ic.v.slice());
    state.base = { m: ic.m.slice(), x: ic.x.slice(), v: ic.v.slice(), lambda: ic.lambda, period: ic.period };
    if (p.tol) setTol(p.tol);
    state.speed = p.speed;
    state.trailDur = p.trail;
    setView3d(!!p.tilt, p.tilt || 0);
    setOpt("rotating", !!p.rotating);
    syncSpeedSlider(); syncTrailSlider();
    $("preset").value = p.id;
    restart(true);
    renderPresetInfo();
    buildMassRows();
    buildICTable();
  }

  function setTol(t) { state.tol = t; $("tol").value = String(t); }

  function restart(fromPreset) {
    const b = state.base;
    const m = b.m, x = b.x.slice(), v = b.v.slice();
    state.sim = new NBody(m, x, v, { tol: state.tol, eps: state.eps });
    state.y0 = Float64Array.from([...x, ...v]);
    state.trails = [0, 1, 2].map(() => new Trail(14000));
    state.twinTrails = [0, 1, 2].map(() => new Trail(6000));
    state.diag = newDiag();
    state.sphere = [];
    state.periods = [];
    state.tDisp = 0;
    state.collision = null;
    state.status = null;
    state.escape = null;
    makeTwin();
    // system size for framing and twin saturation
    let R = 0;
    for (let i = 0; i < 3; i++) R = Math.max(R, Math.hypot(x[3 * i], x[3 * i + 1], x[3 * i + 2]));
    state.sizeScale = Math.max(R, 1e-3);
    state.camR = R * 1.35;
    state.camC = [0, 0, 0];
    state.manualZoom = false;
    state.camPrev = null;
    computeSphereLandmarks();
    const p = state.preset;
    state.period = state.custom ? null : p.massMode === "regen" ? b.period || null : p.period && sameMasses(p.masses) ? p.period : null;
    state.vRef = null; state.aRef = null; state.Lstart = undefined; state.lastZ = undefined; state.crossings = 0;
    if (p.preroll && state.period && !state.custom) { advanceTo(state.period * 1.0001, Infinity); state.tDisp = state.sim.t; }
    sampleDisplay();
    pushDiag(true);
    updateHud(true);
    updateLegend();
  }
  function sameMasses(pm) { return pm && state.base.m.every((q, i) => Math.abs(q - pm[i]) < 1e-12); }

  function makeTwin() {
    state.twin = null;
    if (!state.opts.twin) return;
    const b = state.base, x = b.x.slice(), v = b.v.slice();
    x[0] += state.delta; // body 1 displaced in x
    state.twin = new NBody(b.m, x, v, { tol: state.tol, eps: state.eps });
    state.twinBuf = new Float64Array(18);
    state.diag.S = new Series(1);
    state.diag.lyap = { n: 0, st: 0, sy: 0, stt: 0, sty: 0, horizon: null };
  }

  const tmpX = new Float64Array(9), tmpV = new Float64Array(9);
  function frameAngle(x) { return Math.atan2(x[4] - x[1], x[3] - x[0]); }

  function recordStep(sim, trails) {
    const t0 = sim.prev[18], t1 = sim.t, x0 = sim.prev, x1 = sim.y;
    let dmax = 0;
    for (let i = 0; i < 9; i++) dmax = Math.max(dmax, Math.abs(x1[i] - x0[i]));
    const seg = viewRadius() / 260;
    const n = clamp(Math.ceil(dmax / seg), 1, 48);
    for (let k = 1; k <= n; k++) {
      const tau = t0 + ((t1 - t0) * k) / n;
      if (k === n) for (let i = 0; i < 9; i++) tmpX[i] = x1[i]; else sim.sample(tau, tmpX);
      const th = frameAngle(tmpX);
      const c = Math.cos(th), s = Math.sin(th);
      for (let i = 0; i < 3; i++) trails[i].push(tmpX[3 * i], tmpX[3 * i + 1], tmpX[3 * i + 2], tau, c, s);
    }
  }

  function advanceTo(target, budgetMs) {
    const sim = state.sim, start = performance.now();
    let steps = 0;
    while (sim.t < target) {
      if (state.collision) break;
      try { sim.step(); } catch (e) { state.collision = { t: sim.t }; break; }
      steps++;
      recordStep(sim, state.trails, false);
      afterMainStep();
      if ((steps & 7) === 0 && performance.now() - start > budgetMs) break;
    }
    if (state.twin && !state.collision) {
      const tw = state.twin;
      while (tw.t < Math.min(target, sim.t)) {
        try { tw.step(); } catch (e) { state.twin = null; break; }
        recordStep(tw, state.twinTrails, true);
        if (performance.now() - start > budgetMs * 1.6) break;
      }
    }
    return sim.t >= target;
  }

  function afterMainStep() {
    const sim = state.sim, d = state.diag;
    // energy error at exact step ends (interpolated states would add ~1e-4 noise)
    if (sim.t - d.lastT >= Math.max(1e-3, state.speed / 40) || d.E.t.length < 2) {
      d.lastT = sim.t;
      const dE = sim.relEnergyError();
      d.E.push(sim.t, [Math.log10(Math.max(dE, 1e-17))]);
      const x = sim.y, nn = [0, 1, 2].map((i) => {
        let best = Infinity;
        for (let j = 0; j < 3; j++) if (j !== i) best = Math.min(best, Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]));
        return Math.log10(best);
      });
      d.R.push(sim.t, nn);
    }
    const x = sim.y;
    if (state.preset.id === "sitnikov") {
      const m = state.base.m, z = x[8] - (m[0] * x[2] + m[1] * x[5]) / (m[0] + m[1]);
      if (state.lastZ !== undefined && Math.sign(z) !== Math.sign(state.lastZ) && z !== 0) state.crossings = (state.crossings || 0) + 1;
      state.lastZ = z;
    }
    for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++)
      d.minR = Math.min(d.minR, Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]));
    // twin separation at exactly the same time
    const tw = state.twin;
    if (tw) {
      try {
        while (tw.t < sim.t) { tw.step(); recordStep(tw, state.twinTrails, true); }
        tw.stateAt(sim.t, state.twinBuf);
        let s = 0;
        for (let i = 0; i < 18; i++) { const q = state.twinBuf[i] - sim.y[i]; s += q * q; }
        const sep = Math.sqrt(s);
        const S = d.S;
        if (!S.t.length || sim.t - S.t[S.t.length - 1] >= Math.max(1e-3, state.speed / 40)) {
          S.push(sim.t, [Math.log10(Math.max(sep, 1e-17))]);
          const L = d.lyap, R = state.sizeScale;
          if (sep > 10 * state.delta && sep < 1e-2 * R && sim.t > 0) {
            const y = Math.log(sep); L.n++; L.st += sim.t; L.sy += y; L.stt += sim.t * sim.t; L.sty += sim.t * y;
          }
          if (L.horizon === null && sep > 0.1 * R) L.horizon = sim.t;
        }
      } catch (e) { state.twin = null; }
    }
    // periodic orbit: exact return check at multiples of the period
    if (state.period) {
      const k = Math.floor(sim.t / state.period + 1e-9);
      const kPrev = Math.floor(sim.prev[18] / state.period + 1e-9);
      if (k > kPrev && k >= 1 && k <= 999) {
        const buf = new Float64Array(18);
        sim.stateAt(k * state.period, buf);
        let err = 0;
        for (let i = 0; i < 18; i++) err = Math.max(err, Math.abs(buf[i] - state.y0[i]));
        state.periods.push({ k, err });
      }
    }
  }

  function sampleDisplay() {
    const sim = state.sim;
    const tau = Math.min(state.tDisp, sim.t);
    if (sim.steps === 0) { for (let i = 0; i < 9; i++) { state.disp[i] = sim.y[i]; state.dispV[i] = sim.y[9 + i]; } }
    else sim.sample(Math.max(tau, sim.prev[18]), state.disp, state.dispV);
    state.theta = frameAngle(state.disp);
    const r = [state.disp[3] - state.disp[0], state.disp[4] - state.disp[1]], w = [state.dispV[3] - state.dispV[0], state.dispV[4] - state.dispV[1]];
    state.omega = (r[0] * w[1] - r[1] * w[0]) / Math.max(1e-12, r[0] * r[0] + r[1] * r[1]);
    if (state.twin) {
      const tw = state.twin;
      if (tw.steps === 0) for (let i = 0; i < 9; i++) state.twinDisp[i] = tw.y[i];
      else tw.sample(Math.max(Math.min(tau, tw.t), tw.prev[18]), state.twinDisp);
    }
  }

  function pushDiag() {
    const s = shapeSphere(state.disp, state.base.m);
    state.sphere.push(s);
    if (state.sphere.length > 2500) state.sphere.splice(0, state.sphere.length - 2500);
  }

  /* ================================================================ camera */
  const stage = $("stage"), sky = $("sky"), starsCv = $("stars"), fieldCv = $("field");
  const ctx = sky.getContext("2d");
  let W = 0, H = 0, DPR = 1;
  function viewRadius() { return state.camR || 1; }
  function scale() { return (Math.min(W, H) / 2) / viewRadius(); }
  function toDisplay(x, y, z, c, s) { // inertial -> displayed frame (co-rotating if on)
    if (!state.opts.rotating) return [x, y, z];
    return [c * x + s * y, -s * x + c * y, z];
  }
  function project(x, y, z) {
    const sc = scale(), C = state.camC;
    const dx = x - C[0], dy = y - C[1], dz = z - C[2];
    if (!state.view3d) return [W / 2 + dx * sc, H / 2 - dy * sc, 1, 0];
    const cy = Math.cos(state.yaw), sy = Math.sin(state.yaw);
    const x1 = dx * cy - dy * sy, y1 = dx * sy + dy * cy;
    const T = (state.tilt * Math.PI) / 180, ct = Math.cos(T), st = Math.sin(T);
    const y2 = y1 * ct + dz * st, z2 = -y1 * st + dz * ct;
    const D = 4.2 * viewRadius(), f = D / Math.max(D * 0.15, D - z2);
    return [W / 2 + x1 * sc * f, H / 2 - y2 * sc * f, f, z2];
  }
  function unproject(sx, sy) {
    const sc = scale();
    return [state.camC[0] + (sx - W / 2) / sc, state.camC[1] - (sy - H / 2) / sc];
  }

  function updateCamera(dtReal) {
    if (!state.opts.autofit || state.manualZoom) return;
    const x = state.disp, m = state.base.m;
    const th = state.theta, c = Math.cos(th), s = Math.sin(th);
    const P = [0, 1, 2].map((i) => toDisplay(x[3 * i], x[3 * i + 1], x[3 * i + 2], c, s));
    let idx = [0, 1, 2];
    const st = state.status;
    if (st && st.state === "escape" && st.R > 8 * Math.max(st.a, 0.05)) idx = st.pair;
    let M = 0, C = [0, 0, 0];
    for (const i of idx) { M += m[i]; for (let k = 0; k < 3; k++) C[k] += m[i] * P[i][k]; }
    for (let k = 0; k < 3; k++) C[k] /= M;
    let R = 0;
    for (const i of idx) R = Math.max(R, Math.hypot(P[i][0] - C[0], P[i][1] - C[1], P[i][2] - C[2]));
    if (idx.length === 2 && st) R = Math.max(R, st.a * (1 + st.e) * 0.6);
    // include recent trail extent so loops stay in view
    const tr = state.trails, tNow = state.tDisp, look = Math.min(state.trailDur, 30);
    if (state.opts.trails && idx.length === 3) {
      for (let i = 0; i < 3; i++) {
        const T = tr[i];
        for (let k = 0, stp = Math.max(1, Math.floor(T.len / 300)); k < T.len; k += stp) {
          const j = (T.head - 1 - k + T.cap) % T.cap;
          if (tNow - T.t[j] > look) break;
          const q = toDisplay(T.x[j], T.y[j], T.z[j], T.c[j], T.s[j]);
          R = Math.max(R, Math.hypot(q[0] - C[0], q[1] - C[1], q[2] - C[2]));
        }
      }
    }
    R = Math.max(R * 1.18, state.sizeScale * 0.12);
    const key = idx.join("");
    if (state.camPrev && state.camPrev.key === key) for (let q = 0; q < 3; q++) state.camC[q] += C[q] - state.camPrev.C[q];
    state.camPrev = { key, C };
    const k = 1 - Math.exp(-dtReal * 2.2);
    state.camR += (R - state.camR) * k;
    for (let q = 0; q < 3; q++) state.camC[q] += (C[q] - state.camC[q]) * k;
  }

  /* ================================================================ star field */
  function drawStars() {
    const c = starsCv.getContext("2d");
    starsCv.width = W * DPR; starsCv.height = H * DPR;
    c.setTransform(DPR, 0, 0, DPR, 0, 0);
    c.clearRect(0, 0, W, H);
    const r = mulberry32(42), n = Math.floor((W * H) / 1800);
    for (let i = 0; i < n; i++) {
      const x = r() * W, y = r() * H, b = Math.pow(r(), 3);
      const tint = r();
      c.fillStyle = tint < 0.15 ? `rgba(255,214,170,${0.15 + 0.6 * b})` : tint < 0.3 ? `rgba(170,200,255,${0.15 + 0.6 * b})` : `rgba(225,230,245,${0.12 + 0.55 * b})`;
      const s = 0.4 + 1.1 * b;
      c.beginPath(); c.arc(x, y, s, 0, 2 * Math.PI); c.fill();
    }
  }

  /* ================================================================ gravity field (WebGL) */
  const field = (function () {
    let gl = null;
    try { gl = fieldCv.getContext("webgl", { premultipliedAlpha: true, antialias: false }); } catch (e) { gl = null; }
    if (!gl) return null;
    const deriv = gl.getExtension("OES_standard_derivatives");
    const vs = "attribute vec2 p; void main(){ gl_Position = vec4(p, 0.0, 1.0); }";
    const fs = `${deriv ? "#extension GL_OES_standard_derivatives : enable\n" : ""}precision highp float;
uniform vec2 res; uniform vec2 cen; uniform float sc; uniform vec3 bp[3]; uniform float bm[3];
uniform float soft; uniform float w2; uniform float ref;
void main(){
  vec2 w = cen + (gl_FragCoord.xy - 0.5*res)/sc;
  float phi = 0.0;
  for (int i = 0; i < 3; i++) { vec3 d = vec3(w, 0.0) - bp[i]; phi += bm[i] / sqrt(dot(d,d) + soft*soft); }
  phi += 0.5 * w2 * dot(w, w);
  float v = log(phi / ref);
  float depth = clamp(0.5 + 0.22 * v, 0.0, 1.0);
  vec3 base = mix(vec3(0.02, 0.03, 0.08), vec3(0.20, 0.22, 0.52), depth * depth);
  float a = 0.10 + 0.62 * depth * depth;
  float k = v * 2.6;
  float dist = 0.5 - abs(fract(k) - 0.5);
  ${deriv ? "float fw = max(fwidth(k), 1e-4);" : "float fw = 0.04;"}
  float line = 1.0 - smoothstep(0.0, fw * 1.25, dist);
  line *= smoothstep(-3.5, -1.0, v);
  vec3 col = base + line * vec3(0.42, 0.47, 0.85) * 0.55;
  a = max(a, line * 0.5);
  gl_FragColor = vec4(col * a, a);
}`;
    const sh = (type, src) => { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); return s; };
    let prog;
    try {
      prog = gl.createProgram();
      gl.attachShader(prog, sh(gl.VERTEX_SHADER, vs)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, fs));
      gl.linkProgram(prog);
      if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error("link");
    } catch (e) { return null; }
    gl.useProgram(prog);
    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(prog, "p");
    gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
    const U = {}; ["res", "cen", "sc", "bp", "bm", "soft", "w2", "ref"].forEach((n) => (U[n] = gl.getUniformLocation(prog, n)));
    return {
      draw(P, m, w2) {
        fieldCv.width = Math.round(W * DPR); fieldCv.height = Math.round(H * DPR);
        gl.viewport(0, 0, fieldCv.width, fieldCv.height);
        gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
        gl.uniform2f(U.res, fieldCv.width, fieldCv.height);
        gl.uniform2f(U.cen, state.camC[0], state.camC[1]);
        gl.uniform1f(U.sc, scale() * DPR);
        gl.uniform3fv(U.bp, new Float32Array([...P[0], ...P[1], ...P[2]]));
        gl.uniform1fv(U.bm, new Float32Array(m));
        gl.uniform1f(U.soft, Math.max(state.eps, viewRadius() * 0.004));
        gl.uniform1f(U.w2, w2);
        gl.uniform1f(U.ref, (m[0] + m[1] + m[2]) / viewRadius());
        gl.drawArrays(gl.TRIANGLES, 0, 6);
      },
      clear() { gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT); },
    };
  })();

  /* ================================================================ rendering */
  function resize() {
    const r = stage.getBoundingClientRect();
    W = Math.max(1, Math.round(r.width)); H = Math.max(1, Math.round(r.height));
    DPR = Math.min(window.devicePixelRatio || 1, 2);
    sky.width = Math.round(W * DPR); sky.height = Math.round(H * DPR);
    drawStars();
  }

  function rgba(i, a) { const c = GLOW_RGB[i]; return `rgba(${c[0]},${c[1]},${c[2]},${a})`; }
  function hexA(hex, a) { const n = parseInt(hex.slice(1), 16); return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${a})`; }

  function drawTrails(trails, colorOf, dashed) {
    const B = 14, tNow = state.tDisp, dur = state.trailDur;
    for (let i = 0; i < 3; i++) {
      const T = trails[i];
      if (!T.len) continue;
      const paths = [];
      for (let b = 0; b < B; b++) paths.push(new Path2D());
      let prev = null, prevB = -1;
      const n = T.len, startIdx = (T.head - n + T.cap) % T.cap;
      for (let k = 0; k < n; k++) {
        const j = (startIdx + k) % T.cap;
        const age = tNow - T.t[j];
        if (age > dur) { prev = null; continue; }
        if (age < -1e-9) break;
        const q = toDisplay(T.x[j], T.y[j], T.z[j], T.c[j], T.s[j]);
        const p = project(q[0], q[1], q[2]);
        const b = Math.min(B - 1, Math.floor((1 - age / dur) * B));
        if (prev) {
          if (b !== prevB) { paths[b].moveTo(prev[0], prev[1]); }
          paths[b].lineTo(p[0], p[1]);
        } else paths[b].moveTo(p[0], p[1]);
        prev = p; prevB = b;
      }
      // connect to the live position
      if (prev) {
        const q = dashed ? state.twinDisp : state.disp;
        const c = Math.cos(state.theta), s = Math.sin(state.theta);
        const d = toDisplay(q[3 * i], q[3 * i + 1], q[3 * i + 2], c, s), p = project(d[0], d[1], d[2]);
        paths[B - 1].moveTo(prev[0], prev[1]); paths[B - 1].lineTo(p[0], p[1]);
      }
      ctx.lineCap = "round"; ctx.lineJoin = "round";
      if (dashed) ctx.setLineDash([3, 4]);
      for (let b = 0; b < B; b++) {
        const f = (b + 1) / B;
        ctx.strokeStyle = colorOf(i, dashed ? 0.55 * f : Math.pow(f, 1.35) * 0.95);
        ctx.lineWidth = dashed ? 1 : 0.8 + 1.4 * f;
        ctx.stroke(paths[b]);
      }
      ctx.setLineDash([]);
    }
  }

  function drawGrid(P, m) {
    // rubber-sheet gravity well in the orbital plane (3D view)
    const R = viewRadius() * 1.6, n = 44, ref = (m[0] + m[1] + m[2]) / viewRadius();
    const soft = viewRadius() * 0.08, depth = viewRadius() * 0.24, w2 = state.opts.rotating ? state.omega * state.omega : 0;
    const zAt = (x, y) => {
      let phi = 0;
      for (let i = 0; i < 3; i++) phi += m[i] / Math.sqrt((x - P[i][0]) ** 2 + (y - P[i][1]) ** 2 + P[i][2] ** 2 + soft * soft);
      phi += 0.5 * w2 * (x * x + y * y);
      return -depth * Math.log1p(phi / ref);
    };
    const cx = state.camC[0], cy = state.camC[1];
    ctx.lineWidth = 1;
    for (let pass = 0; pass < 2; pass++) {
      for (let a = 0; a <= n; a++) {
        ctx.beginPath();
        for (let b = 0; b <= n; b++) {
          const u = -R + (2 * R * a) / n, w = -R + (2 * R * b) / n;
          const x = cx + (pass ? w : u), y = cy + (pass ? u : w);
          const p = project(x, y, zAt(x, y));
          if (b === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
        }
        const edge = Math.abs(a - n / 2) / (n / 2);
        ctx.strokeStyle = `rgba(150,165,240,${0.16 * (1 - edge * edge) + 0.03})`;
        ctx.stroke();
      }
    }
  }

  function arrow(x0, y0, x1, y1, color, w) {
    const dx = x1 - x0, dy = y1 - y0, L = Math.hypot(dx, dy);
    if (L < 2) return;
    const ux = dx / L, uy = dy / L, h = Math.min(8, L * 0.35);
    ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = w;
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1 - ux * h * 0.6, y1 - uy * h * 0.6); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - ux * h - uy * h * 0.5, y1 - uy * h + ux * h * 0.5); ctx.lineTo(x1 - ux * h + uy * h * 0.5, y1 - uy * h - ux * h * 0.5); ctx.closePath(); ctx.fill();
  }

  function bodyRadius(i) {
    const m = state.base.m, mmax = Math.max(m[0], m[1], m[2]);
    return clamp(3.2 + 7.5 * Math.cbrt(m[i] / mmax), 2.6, 11);
  }

  function velScale() { // px per velocity unit, so a typical speed draws ~55 px
    let vmax = 1e-9;
    for (let i = 0; i < 3; i++) vmax = Math.max(vmax, Math.hypot(state.dispV[3 * i], state.dispV[3 * i + 1], state.dispV[3 * i + 2]));
    if (!state.vRef || state.editing === null) state.vRef = state.vRef ? Math.max(vmax, state.vRef * 0.98) : vmax;
    return 55 / state.vRef;
  }

  function render() {
    ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
    ctx.clearRect(0, 0, W, H);
    const m = state.base.m, x = state.disp, th = state.theta, c = Math.cos(th), s = Math.sin(th);
    const P = [0, 1, 2].map((i) => toDisplay(x[3 * i], x[3 * i + 1], x[3 * i + 2], c, s));
    // field
    if (field) {
      if (state.opts.field && !state.view3d) field.draw(P, m, state.opts.rotating ? state.omega * state.omega : 0);
      else field.clear();
    }
    if (state.opts.field && state.view3d) drawGrid(P, m);
    // barycentre marker
    const bc = project(0, 0, 0);
    ctx.strokeStyle = "rgba(148,157,181,0.45)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(bc[0] - 5, bc[1]); ctx.lineTo(bc[0] + 5, bc[1]); ctx.moveTo(bc[0], bc[1] - 5); ctx.lineTo(bc[0], bc[1] + 5); ctx.stroke();
    // trails
    if (state.opts.trails) {
      drawTrails(state.trails, (i, a) => hexA(BODY[i], a), false);
      if (state.twin) drawTrails(state.twinTrails, (i, a) => rgba(i, a), true);
    }
    // Trojan / co-rotating: mark the Lagrange points of the 1-2 pair
    if (state.opts.rotating) drawLagrangePoints(P, m);
    // order bodies by depth in 3D
    const order = [0, 1, 2].map((i) => ({ i, p: project(P[i][0], P[i][1], P[i][2]) })).sort((a, b) => a.p[3] - b.p[3]);
    // twin ghosts
    if (state.twin) {
      const tx = state.twinDisp;
      for (let i = 0; i < 3; i++) {
        const q = toDisplay(tx[3 * i], tx[3 * i + 1], tx[3 * i + 2], c, s), p = project(q[0], q[1], q[2]);
        ctx.strokeStyle = rgba(i, 0.85); ctx.lineWidth = 1.3; ctx.setLineDash([2, 2]);
        ctx.beginPath(); ctx.arc(p[0], p[1], bodyRadius(i) * p[2] + 4, 0, 2 * Math.PI); ctx.stroke(); ctx.setLineDash([]);
      }
    }
    const paused = !state.playing;
    const vs = velScale();
    for (const { i, p } of order) {
      const r = bodyRadius(i) * p[2];
      // glow
      ctx.globalCompositeOperation = "lighter";
      const g = ctx.createRadialGradient(p[0], p[1], 0, p[0], p[1], r * 5.5);
      g.addColorStop(0, rgba(i, 0.55)); g.addColorStop(0.25, rgba(i, 0.16)); g.addColorStop(1, rgba(i, 0));
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(p[0], p[1], r * 5.5, 0, 2 * Math.PI); ctx.fill();
      ctx.globalCompositeOperation = "source-over";
      const g2 = ctx.createRadialGradient(p[0] - r * 0.3, p[1] - r * 0.3, 0, p[0], p[1], r);
      g2.addColorStop(0, "#ffffff"); g2.addColorStop(0.45, GLOW[i]); g2.addColorStop(1, BODY[i]);
      ctx.fillStyle = g2; ctx.beginPath(); ctx.arc(p[0], p[1], r, 0, 2 * Math.PI); ctx.fill();
      ctx.strokeStyle = "rgba(6,9,18,0.9)"; ctx.lineWidth = 1.2; ctx.stroke();
      // vectors
      const v = toDisplay(state.dispV[3 * i], state.dispV[3 * i + 1], state.dispV[3 * i + 2], c, s);
      if (state.opts.velocity || (paused && !state.view3d && !state.opts.rotating)) {
        const tip = velocityTip(i, P, vs);
        arrow(p[0], p[1], tip[0], tip[1], hexA(GLOW[i], 0.9), 1.6);
        if (paused && !state.view3d && !state.opts.rotating) { ctx.strokeStyle = hexA(GLOW[i], 0.9); ctx.lineWidth = 1.2; ctx.beginPath(); ctx.arc(tip[0], tip[1], 5, 0, 2 * Math.PI); ctx.stroke(); }
      }
      if (state.opts.forces) {
        const a = accelOf(i, x, m);
        const ad = toDisplay(a[0], a[1], a[2], c, s), an = Math.hypot(ad[0], ad[1], ad[2]);
        const L = 26 * Math.log1p(an / (state.aRef || (state.aRef = an || 1)) * 3);
        const q = project(P[i][0] + (ad[0] / an) * 1e-3 * viewRadius(), P[i][1] + (ad[1] / an) * 1e-3 * viewRadius(), P[i][2] + (ad[2] / an) * 1e-3 * viewRadius());
        const ux = q[0] - p[0], uy = q[1] - p[1], ul = Math.hypot(ux, uy) || 1;
        arrow(p[0], p[1], p[0] + (ux / ul) * L, p[1] + (uy / ul) * L, "rgba(230,233,242,0.75)", 1.2);
      }
      if (state.opts.labels) {
        ctx.font = "500 11px 'IBM Plex Mono', ui-monospace, monospace";
        ctx.fillStyle = "rgba(230,233,242,0.92)";
        ctx.fillText(String(i + 1), p[0] + r + 5, p[1] - r - 2);
      }
    }
    drawOffscreen(P);
    if (state.collision) {
      ctx.font = "500 13px 'IBM Plex Sans', system-ui, sans-serif"; ctx.fillStyle = "#ef7468";
      ctx.fillText(`Solver stopped at t = ${state.collision.t.toFixed(3)}: an exact collision. Add softening ε to pass through.`, 16, H - 46);
    }
  }

  function velocityTip(i, P, vs) {
    const c = Math.cos(state.theta), s = Math.sin(state.theta);
    const v = toDisplay(state.dispV[3 * i], state.dispV[3 * i + 1], state.dispV[3 * i + 2], c, s);
    const p = project(P[i][0], P[i][1], P[i][2]);
    const q = project(P[i][0] + v[0] * 1e-3, P[i][1] + v[1] * 1e-3, P[i][2] + v[2] * 1e-3);
    const vx = (q[0] - p[0]) / 1e-3, vy = (q[1] - p[1]) / 1e-3, sc = scale();
    return [p[0] + (vx / sc) * vs, p[1] + (vy / sc) * vs];
  }

  function accelOf(i, x, m) {
    const a = [0, 0, 0], e2 = state.eps * state.eps;
    for (let j = 0; j < 3; j++) {
      if (j === i) continue;
      const d = [x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]];
      const r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + e2, f = m[j] / (r2 * Math.sqrt(r2));
      for (let k = 0; k < 3; k++) a[k] += f * d[k];
    }
    return a;
  }

  function drawLagrangePoints(P, m) {
    // L1–L5 of the (1, 2) pair, treating body 3 as a test mass: zeros of the effective force on the x-axis, and the two triangles
    const m1 = m[0], m2 = m[1], M = m1 + m2;
    const d = Math.hypot(P[1][0] - P[0][0], P[1][1] - P[0][1]);
    if (!(d > 0)) return;
    const ux = (P[1][0] - P[0][0]) / d, uy = (P[1][1] - P[0][1]) / d;
    const cx = (m1 * P[0][0] + m2 * P[1][0]) / M, cy = (m1 * P[0][1] + m2 * P[1][1]) / M;
    const x1 = -m2 / M * d, x2 = m1 / M * d, w2 = M / (d * d * d);
    const F = (x) => -m1 * (x - x1) / Math.abs(x - x1) ** 3 - m2 * (x - x2) / Math.abs(x - x2) ** 3 + w2 * x;
    const root = (a, b) => { let fa = F(a); for (let k = 0; k < 80; k++) { const mid = 0.5 * (a + b), fm = F(mid); if (Math.sign(fm) === Math.sign(fa)) { a = mid; fa = fm; } else b = mid; } return 0.5 * (a + b); };
    const L = [
      ["L1", root(x1 + 1e-6 * d, x2 - 1e-6 * d), 0], ["L2", root(x2 + 1e-6 * d, x2 + 2 * d), 0], ["L3", root(x1 - 2 * d, x1 - 1e-6 * d), 0],
      ["L4", (x1 + x2) / 2, (Math.sqrt(3) / 2) * d], ["L5", (x1 + x2) / 2, (-Math.sqrt(3) / 2) * d],
    ];
    ctx.font = "500 10px 'IBM Plex Mono', ui-monospace, monospace";
    for (const [name, a, b] of L) {
      const X = cx + a * ux - b * uy, Y = cy + a * uy + b * ux;
      const p = project(X, Y, 0);
      ctx.strokeStyle = "rgba(170,184,255,0.75)"; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(p[0] - 4, p[1] - 4); ctx.lineTo(p[0] + 4, p[1] + 4); ctx.moveTo(p[0] + 4, p[1] - 4); ctx.lineTo(p[0] - 4, p[1] + 4); ctx.stroke();
      ctx.fillStyle = "rgba(170,184,255,0.9)"; ctx.fillText(name, p[0] + 6, p[1] - 6);
    }
  }

  function drawOffscreen(P) {
    const pad = 18;
    for (let i = 0; i < 3; i++) {
      const p = project(P[i][0], P[i][1], P[i][2]);
      if (p[0] > -4 && p[0] < W + 4 && p[1] > -4 && p[1] < H + 4) continue;
      const cx = W / 2, cy = H / 2, dx = p[0] - cx, dy = p[1] - cy;
      const k = Math.min((W / 2 - pad) / Math.abs(dx || 1e-9), (H / 2 - pad) / Math.abs(dy || 1e-9));
      const ex = cx + dx * k, ey = cy + dy * k, ang = Math.atan2(dy, dx);
      ctx.save(); ctx.translate(ex, ey); ctx.rotate(ang);
      ctx.fillStyle = GLOW[i]; ctx.beginPath(); ctx.moveTo(8, 0); ctx.lineTo(-5, -6); ctx.lineTo(-5, 6); ctx.closePath(); ctx.fill();
      ctx.restore();
      const dist = Math.hypot(P[i][0] - state.camC[0], P[i][1] - state.camC[1], P[i][2] - state.camC[2]);
      ctx.font = "500 11px 'IBM Plex Mono', ui-monospace, monospace"; ctx.fillStyle = "rgba(230,233,242,0.9)";
      const label = `${i + 1} · r ${fmt(dist, 1)}`;
      const tw = ctx.measureText(label).width;
      ctx.fillText(label, clamp(ex - tw / 2 - Math.cos(ang) * 30, 6, W - tw - 6), clamp(ey - Math.sin(ang) * 22 + 4, 14, H - 8));
    }
  }

  /* ================================================================ charts */
  const tooltip = $("tooltip");
  function showTip(evt, html) { tooltip.innerHTML = html; tooltip.hidden = false; tooltip.style.left = evt.clientX + "px"; tooltip.style.top = evt.clientY + "px"; }
  function hideTip() { tooltip.hidden = true; }

  class MiniChart {
    constructor(canvas, opts) {
      this.cv = canvas; this.opts = opts; this.hover = null;
      canvas.addEventListener("pointermove", (e) => { const r = canvas.getBoundingClientRect(); this.hover = (e.clientX - r.left) / r.width; this.evt = e; });
      canvas.addEventListener("pointerleave", () => { this.hover = null; hideTip(); });
    }
    draw(series) {
      const cv = this.cv, r = cv.getBoundingClientRect(), w = Math.max(10, r.width), h = Math.max(10, r.height);
      if (cv.width !== Math.round(w * DPR)) { cv.width = Math.round(w * DPR); cv.height = Math.round(h * DPR); }
      const c = cv.getContext("2d");
      c.setTransform(DPR, 0, 0, DPR, 0, 0); c.clearRect(0, 0, w, h);
      const L = 34, Rm = 8, T = 8, B = 16, pw = w - L - Rm, ph = h - T - B;
      if (!series || !series.t.length) {
        c.fillStyle = "#5d6684"; c.font = "400 11px 'IBM Plex Sans', system-ui, sans-serif";
        c.fillText(this.opts.empty || "no data yet", L, T + ph / 2 + 4);
        return;
      }
      const t = series.t, tmax = Math.max(t[t.length - 1], 1e-9), tmin = Math.min(t[0], 0);
      let lo = Infinity, hi = -Infinity;
      for (const ys of series.y) for (const y of ys) { if (isFinite(y)) { lo = Math.min(lo, y); hi = Math.max(hi, y); } }
      if (this.opts.yMin !== undefined) lo = Math.max(lo, this.opts.yMin);
      lo = Math.floor(lo - 0.3); hi = Math.ceil(hi + 0.3);
      if (hi - lo < 2) { hi += 1; lo -= 1; }
      const X = (tt) => L + ((tt - tmin) / (tmax - tmin)) * pw, Y = (yy) => T + (1 - (yy - lo) / (hi - lo)) * ph;
      // grid + ticks (log10 values)
      c.font = "400 9.5px 'IBM Plex Mono', ui-monospace, monospace"; c.textAlign = "right";
      const step = Math.max(1, Math.ceil((hi - lo) / 4));
      for (let yy = Math.ceil(lo / step) * step; yy <= hi; yy += step) {
        c.strokeStyle = "rgba(42,54,87,0.55)"; c.lineWidth = 1; c.beginPath(); c.moveTo(L, Y(yy) + 0.5); c.lineTo(L + pw, Y(yy) + 0.5); c.stroke();
        c.fillStyle = "#7c86a2"; c.fillText(`1e${yy}`, L - 4, Y(yy) + 3);
      }
      c.textAlign = "left"; c.fillStyle = "#7c86a2";
      c.fillText(`t = ${fmt(tmin, 0)}`, L, h - 4);
      c.textAlign = "right"; c.fillText(`${fmt(tmax, 1)}`, L + pw, h - 4); c.textAlign = "left";
      // lines
      series.y.forEach((ys, k) => {
        c.strokeStyle = this.opts.colors[k]; c.lineWidth = this.opts.lineWidth || 1.6; c.lineJoin = "round";
        c.beginPath();
        let started = false;
        for (let j = 0; j < t.length; j++) {
          if (!isFinite(ys[j])) { started = false; continue; }
          const px = X(t[j]), py = Y(clamp(ys[j], lo, hi));
          if (!started) { c.moveTo(px, py); started = true; } else c.lineTo(px, py);
        }
        c.stroke();
        const j = t.length - 1;
        if (isFinite(ys[j])) { c.fillStyle = this.opts.colors[k]; c.beginPath(); c.arc(X(t[j]), Y(clamp(ys[j], lo, hi)), 2.5, 0, 2 * Math.PI); c.fill(); }
      });
      if (this.opts.overlay) this.opts.overlay(c, X, Y, lo, hi);
      // hover crosshair
      if (this.hover !== null) {
        const tt = tmin + clamp((this.hover * w - L) / pw, 0, 1) * (tmax - tmin);
        let j = 0, lo2 = 0, hi2 = t.length - 1;
        while (lo2 < hi2) { const mid = (lo2 + hi2) >> 1; if (t[mid] < tt) lo2 = mid + 1; else hi2 = mid; }
        j = lo2;
        c.strokeStyle = "rgba(230,233,242,0.5)"; c.beginPath(); c.moveTo(X(t[j]) + 0.5, T); c.lineTo(X(t[j]) + 0.5, T + ph); c.stroke();
        const rows = series.y.map((ys, k) => `<div><span style="color:${this.opts.colors[k]}">●</span> <span class="k">${this.opts.names[k]}</span> ${isFinite(ys[j]) ? sci(Math.pow(10, ys[j])) : "–"}</div>`).join("");
        showTip(this.evt, `<div class="k">t = ${fmt(t[j], 3)}</div>${rows}`);
      }
    }
  }
  const chartE = new MiniChart($("chartE"), { colors: ["#e6e9f2"], names: ["|ΔE/E₀|"], yMin: -17 });
  const chartR = new MiniChart($("chartR"), { colors: BODY, names: ["body 1", "body 2", "body 3"] });
  const chartS = new MiniChart($("chartS"), {
    colors: ["#aab8ff"], names: ["separation"], empty: "Turn on the butterfly twin (B) to measure chaos",
    overlay: (c, X, Y, lo, hi) => {
      const L = state.diag.lyap;
      if (L.n < 6) return;
      const lam = (L.n * L.sty - L.st * L.sy) / (L.n * L.stt - L.st * L.st), b = (L.sy - lam * L.st) / L.n;
      const S = state.diag.S, t0 = S.t[0], t1 = S.t[S.t.length - 1];
      c.setLineDash([4, 3]); c.strokeStyle = "rgba(148,157,181,0.8)"; c.lineWidth = 1;
      c.beginPath();
      const ya = (lam * t0 + b) / Math.LN10, yb = (lam * t1 + b) / Math.LN10;
      c.moveTo(X(t0), Y(clamp(ya, lo, hi))); c.lineTo(X(t1), Y(clamp(yb, lo, hi))); c.stroke(); c.setLineDash([]);
    },
  });
  $("chartRLegend").innerHTML = [0, 1, 2].map((i) => `<span class="legend-item"><span class="dot" style="background:${BODY[i]}"></span>body ${i + 1}</span>`).join("");

  /* ================================================================ shape sphere */
  const sphereCv = $("sphere");
  function computeSphereLandmarks() {
    const m = state.base.m;
    const coll = [[0, 1], [1, 2], [0, 2]].map(([i, j]) => {
      const x = [0, 0, 0, 0, 0, 0, 0, 0, 0];
      const k = 3 - i - j; x[3 * k] = 1; // the pair (i, j) sits together at the origin
      return { pt: shapeSphere(x, m), name: `${i + 1}${j + 1}` };
    });
    const eul = [0, 1, 2].map((mid) => {
      const order = [0, 1, 2].filter((q) => q !== mid); order.splice(1, 0, mid);
      const mm = order.map((q) => m[q]);
      const e = eulerIC(mm, 2), x = new Array(9).fill(0);
      order.forEach((q, slot) => { x[3 * q] = e.x[3 * slot]; });
      return shapeSphere(x, m);
    });
    state.sphereMarks = { coll, eul };
  }
  function drawSphere() {
    const cv = sphereCv, r0 = cv.getBoundingClientRect(), S = Math.max(40, r0.width);
    if (cv.width !== Math.round(S * DPR)) { cv.width = Math.round(S * DPR); cv.height = Math.round(S * DPR); }
    const c = cv.getContext("2d");
    c.setTransform(DPR, 0, 0, DPR, 0, 0); c.clearRect(0, 0, S, S);
    const R = S / 2 - 6, cx = S / 2, cy = S / 2;
    const tilt = 0.42 + 0.0 * Math.sin(performance.now() / 4000), spin = performance.now() / 9000 * (reduceMotion ? 0 : 1);
    const P = (w) => { // w = (w1, w2, w3); w3 is the pole axis
      const x = w[0] * Math.cos(spin) - w[1] * Math.sin(spin), y = w[0] * Math.sin(spin) + w[1] * Math.cos(spin), z = w[2];
      const y2 = y * Math.sin(tilt) + z * Math.cos(tilt), d = y * Math.cos(tilt) - z * Math.sin(tilt);
      return [cx + x * R, cy - y2 * R, -d];
    };
    const g = c.createRadialGradient(cx - R * 0.3, cy - R * 0.35, R * 0.1, cx, cy, R);
    g.addColorStop(0, "rgba(40,50,92,0.9)"); g.addColorStop(1, "rgba(10,14,30,0.95)");
    c.fillStyle = g; c.beginPath(); c.arc(cx, cy, R, 0, 2 * Math.PI); c.fill();
    c.strokeStyle = "rgba(42,54,87,1)"; c.lineWidth = 1; c.stroke();
    const ring = (f, col, w) => { // latitude/meridian line
      c.strokeStyle = col; c.lineWidth = w;
      let prev = null;
      for (let k = 0; k <= 96; k++) {
        const p = P(f((2 * Math.PI * k) / 96));
        if (prev && prev[2] >= 0 && p[2] >= 0) { c.beginPath(); c.moveTo(prev[0], prev[1]); c.lineTo(p[0], p[1]); c.stroke(); }
        prev = p;
      }
    };
    for (const lat of [-0.5, 0.5]) ring((a) => [Math.cos(a) * Math.cos(lat), Math.sin(a) * Math.cos(lat), Math.sin(lat)], "rgba(60,74,118,0.6)", 0.8);
    for (let mer = 0; mer < 3; mer++) { const ph = (mer * Math.PI) / 3; ring((a) => [Math.cos(a) * Math.cos(ph), Math.cos(a) * Math.sin(ph), Math.sin(a)], "rgba(60,74,118,0.45)", 0.8); }
    ring((a) => [Math.cos(a), Math.sin(a), 0], "rgba(148,157,181,0.9)", 1.2);
    // path
    const pts = state.sphere, n = pts.length;
    for (let k = 1; k < n; k++) {
      const a = P(pts[k - 1]), b = P(pts[k]);
      const front = a[2] >= 0 && b[2] >= 0, f = k / n;
      c.strokeStyle = `rgba(230,233,242,${(front ? 0.85 : 0.22) * (0.15 + 0.85 * f)})`; c.lineWidth = front ? 1.3 : 0.8;
      c.beginPath(); c.moveTo(a[0], a[1]); c.lineTo(b[0], b[1]); c.stroke();
    }
    // landmarks
    const mk = state.sphereMarks;
    for (const pole of [[0, 0, 1], [0, 0, -1]]) { const p = P(pole); c.fillStyle = p[2] >= 0 ? "#aab8ff" : "rgba(170,184,255,0.3)"; c.beginPath(); c.arc(p[0], p[1], 3, 0, 2 * Math.PI); c.fill(); }
    for (const e of mk.eul) { const p = P(e); c.fillStyle = p[2] >= 0 ? "#949db5" : "rgba(148,157,181,0.3)"; c.beginPath(); c.arc(p[0], p[1], 2.4, 0, 2 * Math.PI); c.fill(); }
    for (const q of mk.coll) {
      const p = P(q.pt); c.strokeStyle = p[2] >= 0 ? "#ef7468" : "rgba(239,116,104,0.3)"; c.lineWidth = 1.4;
      c.beginPath(); c.moveTo(p[0] - 3, p[1] - 3); c.lineTo(p[0] + 3, p[1] + 3); c.moveTo(p[0] + 3, p[1] - 3); c.lineTo(p[0] - 3, p[1] + 3); c.stroke();
    }
    if (n) { const p = P(pts[n - 1]); c.fillStyle = "#ffffff"; c.beginPath(); c.arc(p[0], p[1], 3.2, 0, 2 * Math.PI); c.fill(); c.strokeStyle = "rgba(6,9,18,0.9)"; c.stroke(); }
  }

  /* ================================================================ HUD + readouts */
  function updateLegend() {
    const m = state.base.m;
    $("hudLegend").innerHTML = [0, 1, 2].map((i) => `<span class="legend-item"><span class="dot" style="background:${GLOW[i]};box-shadow:0 0 8px ${GLOW[i]}"></span><b>${i + 1}</b> m = ${fmtMass(m[i])}</span>`).join("");
  }
  function fmtMass(m) { return m >= 0.01 ? (+m.toPrecision(3)).toString() : m.toExponential(2); }

  function classifyNow() {
    const p = state.preset;
    const c = classify(state.sim.y, state.sim.y.subarray(9), state.base.m, 1, p.rMin || 0);
    state.status = c.state === "unbound" ? { state: "unbound" } : c;
    if (c.state === "escape" && !state.escape) state.escape = { t: state.sim.t, k: c.third };
  }

  function updateHud(force) {
    const sim = state.sim, p = state.preset;
    $("hudName").textContent = state.custom ? `${p.name} (edited)` : p.name;
    let t = `t = <b>${state.tDisp.toFixed(2)}</b>`;
    if (p.units === "days" && !state.custom) t += ` · ${(state.tDisp * DAY_PER_TU).toFixed(0)} days`;
    if (state.period) t += ` · period ${(state.tDisp / state.period).toFixed(2)}`;
    $("hudTime").innerHTML = t;
    // status pill
    const pill = $("pill"), st = state.status;
    let kind = "interplay", text = "Three-way interplay";
    if (state.collision) { kind = "collision"; text = "Exact collision: singular point"; }
    else if (state.period && (!state.periods.length || state.periods[state.periods.length - 1].err < 1e-3)) {
      kind = "periodic"; const last = state.periods[state.periods.length - 1];
      text = last ? `Periodic · back within ${sci(last.err)} after ${last.k} ${last.k === 1 ? "period" : "periods"}` : `Periodic orbit · T = ${state.period.toFixed(4)}`;
    } else if (st) {
      if (st.state === "escape") { kind = "escape"; text = `Body ${st.third + 1} escaping · v∞ ${fmt(st.vInf, 2)} · ${st.pair[0] + 1}–${st.pair[1] + 1} binary, e ${fmt(st.e, 2)}`; }
      else if (st.state === "hierarchical") { kind = "hierarchical"; text = `Hierarchical · ${st.pair[0] + 1}–${st.pair[1] + 1} pair (e ${fmt(st.e, 2)}) with ${st.third + 1} orbiting it`; }
      else if (st.state === "unbound") { kind = "escape"; text = "Unbound · all three flying apart"; }
      if (state.period && state.periods.length) text = `Orbit broke after ${state.periods.length} periods · ` + text.split(" · ")[0];
    }
    if (p.id === "trojan" && !state.custom && !state.collision) {
      const x = state.disp, a2 = Math.atan2(x[4] - x[1], x[3] - x[0]), a3 = Math.atan2(x[7] - x[1], x[6] - x[0]);
      let d = ((a3 - a2) * 180) / Math.PI; d = ((d + 540) % 360) - 180;
      kind = "periodic"; text = `Trojan librating about L4 · ${d.toFixed(0)}° ahead of the planet (L4 = 60°)`;
    }
    if (p.id === "sitnikov" && !state.custom && !state.collision && (!st || st.state !== "escape")) {
      const x = state.disp, m = state.base.m, z = x[8] - (m[0] * x[2] + m[1] * x[5]) / (m[0] + m[1]);
      kind = "interplay"; text = `Body 3 bobbing · z = ${z >= 0 ? "+" : "−"}${Math.abs(z).toFixed(2)} · ${state.crossings} plane crossings`;
    }
    pill.dataset.state = kind; $("pillText").textContent = text;
    $("stateWord").textContent = kind === "periodic" ? "periodic" : st ? st.state : "";
    // readouts
    const E = sim.energy(), L = sim.angularMomentum();
    $("rE").textContent = fmt(E, 6);
    $("rdE").textContent = sci(sim.relEnergyError());
    const Ln = Math.hypot(L[0], L[1], L[2]);
    if (sim.steps === 0 || state.Lstart === undefined || force) state.Lstart = Ln;
    $("rL").textContent = fmt(Ln, 6);
    $("rdL").textContent = sci(Math.abs(Ln - state.Lstart));
    $("rSteps").textContent = sim.steps.toLocaleString();
    $("rEvals").textContent = sim.forceEvals.toLocaleString();
    $("rDt").textContent = sim.steps ? sci(sim.t - sim.prev[18]) : "–";
    $("rMin").textContent = isFinite(state.diag.minR) ? fmt(state.diag.minR, 4) : "–";
    if (st && st.pair) {
      $("rPair").textContent = `${st.pair[0] + 1}–${st.pair[1] + 1}${state.view3d ? ` · i ${st.inc.toFixed(0)}°` : ""}`;
      $("rAE").textContent = `${fmt(st.a, 3)} · ${fmt(st.e, 3)}`;
    } else { $("rPair").textContent = "–"; $("rAE").textContent = "–"; }
    const Ld = state.diag.lyap;
    if (state.twin && Ld.n >= 6) {
      const lam = (Ld.n * Ld.sty - Ld.st * Ld.sy) / (Ld.n * Ld.stt - Ld.st * Ld.st);
      const hz = Ld.horizon !== null ? ` · diverged at t ${Ld.horizon.toFixed(1)}` : "";
      $("lyap").textContent = `λ ≈ ${lam.toFixed(3)}${hz}`;
    } else $("lyap").textContent = state.opts.twin ? `δ = ${sci(state.delta)} · measuring` : "turn on the twin";
    $("periodRead").textContent = state.period && state.periods.length ? `return error ${sci(state.periods[state.periods.length - 1].err)}` : "";
  }

  function renderPresetInfo() {
    const p = state.preset;
    let b = p.blurb;
    if (p.routh) {
      const r = routh(state.base.m);
      b += r.stable ? `<b>Routh: stable</b> (27Σmᵢmⱼ = ${fmt(r.lhs, 2)} < (Σm)² = ${fmt(r.rhs, 2)}).` : `<b>Routh: unstable</b> for these masses (27Σmᵢmⱼ = ${fmt(r.lhs, 2)} ≥ (Σm)² = ${fmt(r.rhs, 2)}): rounding errors will grow until it breaks.`;
    }
    if (p.euler && state.base.lambda) b += ` Here λ = ${state.base.lambda.toFixed(6)}.`;
    $("blurb").innerHTML = b;
    $("source").textContent = p.source ? `Source: ${p.source}` : "";
    $("randomRow").hidden = !p.random;
    $("seedLabel").textContent = p.random ? `seed ${state.seed}` : "";
    $("massNote").textContent = p.scaleNote ? `masses in ${p.scaleNote.split(", ")[1]}, G = 1` : p.massMode === "regen" ? "drag: exact solution rebuilds" : "masses, G = 1";
  }

  /* ================================================================ controls */
  const MASS_MIN = -3, MASS_MAX = 1.3;
  function buildMassRows() {
    const box = $("massRows"); box.innerHTML = "";
    state.base.m.forEach((mv, i) => {
      const row = document.createElement("div"); row.className = "mass";
      const id = `mass${i}`;
      const lv = clamp(Math.log10(Math.max(mv, 1e-9)), MASS_MIN, MASS_MAX);
      row.innerHTML = `<span class="dot" style="background:${GLOW[i]}" aria-hidden="true"></span><input type="range" id="${id}" min="${MASS_MIN}" max="${MASS_MAX}" step="0.01" value="${lv}" aria-label="Mass of body ${i + 1} (log scale)"><output for="${id}" id="${id}o">${fmtMass(mv)}</output>`;
      box.appendChild(row);
      const inp = row.querySelector("input");
      inp.addEventListener("input", () => {
        const val = Math.pow(10, +inp.value);
        state.masses[i] = +val.toPrecision(3);
        row.querySelector("output").textContent = fmtMass(state.masses[i]);
        applyMasses();
      });
    });
  }
  function applyMasses() {
    const p = state.preset, wasPlaying = state.playing;
    if (p.massMode === "regen") {
      const ic = p.make(state.masses.slice());
      state.base = { m: ic.m.slice(), x: ic.x.slice(), v: ic.v.slice(), lambda: ic.lambda, period: ic.period };
    } else {
      const ic = toCom(state.masses.slice(), state.base.x.slice(), state.base.v.slice());
      state.base = { m: ic.m, x: ic.x, v: ic.v };
      state.custom = !sameMasses(p.masses);
    }
    restart();
    state.playing = wasPlaying; syncPlay();
    renderPresetInfo(); buildICTable();
  }

  function buildICTable() {
    const tb = $("icTable"), b = state.base;
    const cols = ["x", "y", "z", "vx", "vy", "vz"];
    let h = `<thead><tr><th></th>${cols.map((c) => `<th>${c}</th>`).join("")}</tr></thead><tbody>`;
    for (let i = 0; i < 3; i++) {
      h += `<tr><th><span class="dot" style="background:${GLOW[i]};display:inline-block"></span></th>`;
      cols.forEach((c, k) => {
        const val = k < 3 ? b.x[3 * i + k] : b.v[3 * i + k - 3];
        h += `<td><input type="number" step="0.01" id="ic${i}${c}" aria-label="Body ${i + 1} ${c}" value="${+val.toPrecision(6)}"></td>`;
      });
      h += "</tr>";
    }
    tb.innerHTML = h + "</tbody>";
    tb.querySelectorAll("input").forEach((inp) => inp.addEventListener("change", () => {
      const x = [], v = [];
      for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) {
        x.push(+$(`ic${i}${cols[k]}`).value || 0); v.push(+$(`ic${i}${cols[k + 3]}`).value || 0);
      }
      setCustom(x, v);
    }));
  }

  function setCustom(x, v) {
    const ic = toCom(state.base.m.slice(), x, v);
    state.base = { m: ic.m, x: ic.x, v: ic.v };
    state.custom = true;
    const wasPlaying = state.playing;
    restart();
    state.playing = wasPlaying; syncPlay();
    buildICTable();
  }

  function populatePresetSelect() {
    const sel = $("preset");
    const groups = {};
    for (const p of PRESETS) (groups[p.group] = groups[p.group] || []).push(p);
    sel.innerHTML = Object.entries(groups).map(([g, ps]) => `<optgroup label="${g}">${ps.map((p) => `<option value="${p.id}">${p.name}</option>`).join("")}</optgroup>`).join("");
    sel.addEventListener("change", () => { loadPreset(sel.value); });
  }

  function speedFromSlider(v) { return Math.pow(10, -1.5 + (3.9 * v) / 100); }
  function syncSpeedSlider() {
    $("speed").value = String(clamp(((Math.log10(state.speed) + 1.5) / 3.9) * 100, 0, 100));
    $("speedOut").textContent = `${state.speed >= 10 ? state.speed.toFixed(0) : state.speed.toFixed(2)} t/s`;
  }
  $("speed").addEventListener("input", () => { state.speed = speedFromSlider(+$("speed").value); $("speedOut").textContent = `${state.speed >= 10 ? state.speed.toFixed(0) : state.speed.toFixed(2)} t/s`; });

  function trailFromSlider(v) { return Math.pow(10, -0.5 + (3 * v) / 100); }
  function syncTrailSlider() { $("trailLen").value = String(clamp(((Math.log10(state.trailDur) + 0.5) / 3) * 100, 0, 100)); $("trailOut").textContent = fmt(state.trailDur, state.trailDur < 10 ? 1 : 0); }
  $("trailLen").addEventListener("input", () => { state.trailDur = trailFromSlider(+$("trailLen").value); $("trailOut").textContent = fmt(state.trailDur, state.trailDur < 10 ? 1 : 0); });

  function setOpt(name, val) {
    state.opts[name] = val;
    const b = document.querySelector(`.chip[data-opt="${name}"]`);
    if (b) b.setAttribute("aria-pressed", String(val));
    if (name === "twin" && state.sim) {
      if (val) { makeTwin(); state.twinTrails.forEach((t) => t.clear());
        // catch the twin up to the present
        if (state.twin && state.sim.t > 0) { const tgt = state.sim.t; while (state.twin && state.twin.t < tgt) { state.twin.step(); recordStep(state.twin, state.twinTrails, true); } state.diag.S = new Series(1); }
      } else { state.twin = null; state.twinTrails.forEach((t) => t.clear()); }
    }
    if (name === "autofit" && val) state.manualZoom = false;
  }
  document.querySelectorAll(".chip").forEach((b) => b.addEventListener("click", () => setOpt(b.dataset.opt, b.getAttribute("aria-pressed") !== "true")));

  function setView3d(on, tilt) {
    state.view3d = on;
    state.tilt = on ? (tilt || 60) : 0;
    if (!on) state.yaw = 0;
    $("viewPlane").setAttribute("aria-pressed", String(!on));
    $("view3d").setAttribute("aria-pressed", String(on));
  }
  $("viewPlane").addEventListener("click", () => setView3d(false));
  $("view3d").addEventListener("click", () => setView3d(true, state.preset.tilt || 60));

  function syncPlay() {
    $("playLabel").textContent = state.playing ? "Pause" : "Play";
    $("playIcon").setAttribute("d", state.playing ? "M3 2h3v10H3zM8 2h3v10H8z" : "M3 2l9 5-9 5z");
    $("hint").hidden = state.playing || state.view3d || state.opts.rotating;
    stage.classList.toggle("editing", !state.playing && !state.view3d);
  }
  $("play").addEventListener("click", () => { state.playing = !state.playing; syncPlay(); });
  $("restart").addEventListener("click", () => { restart(); });
  $("reroll").addEventListener("click", () => { state.seed = (state.seed * 7919 + 13) % 100000; loadPreset("random"); });
  $("tol").addEventListener("change", () => { state.tol = +$("tol").value; restart(); });
  $("eps").addEventListener("change", () => { state.eps = +$("eps").value; restart(); });
  $("delta").addEventListener("change", () => { state.delta = +$("delta").value; if (state.opts.twin) { setOpt("twin", false); setOpt("twin", true); } });

  document.addEventListener("keydown", (e) => {
    if (e.target && /INPUT|SELECT|TEXTAREA/.test(e.target.tagName)) return;
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    const k = e.key.toLowerCase();
    if (k === " " && e.target && e.target.tagName === "BUTTON") return; // the focused button handles Space itself
    if (k === " ") { e.preventDefault(); state.playing = !state.playing; syncPlay(); }
    else if (k === "r") restart();
    else if (k === "b") setOpt("twin", !state.opts.twin);
    else if (k === "f") setOpt("field", !state.opts.field);
    else if (k === "l") setOpt("trails", !state.opts.trails);
    else if (k === "c") setOpt("rotating", !state.opts.rotating);
    else if (k === "v") setView3d(!state.view3d, state.preset.tilt || 60);
  });

  /* ---- pointer: orbit / pan / zoom / edit */
  let drag = null;
  stage.addEventListener("pointerdown", (e) => {
    if (e.target.closest(".hud > *") && e.target !== sky) return;
    const r = stage.getBoundingClientRect(), sx = e.clientX - r.left, sy = e.clientY - r.top;
    stage.setPointerCapture(e.pointerId);
    // edit mode: paused, top view, inertial frame
    if (!state.playing && !state.view3d && !state.opts.rotating) {
      const x = state.disp, P = [0, 1, 2].map((i) => [x[3 * i], x[3 * i + 1], x[3 * i + 2]]), vs = velScale();
      for (let i = 0; i < 3; i++) {
        const tip = velocityTip(i, P, vs);
        if (Math.hypot(sx - tip[0], sy - tip[1]) < 11) { drag = { mode: "vel", i, vs }; state.editing = i; return; }
      }
      for (let i = 0; i < 3; i++) {
        const p = project(P[i][0], P[i][1], P[i][2]);
        if (Math.hypot(sx - p[0], sy - p[1]) < Math.max(12, bodyRadius(i) + 6)) { drag = { mode: "pos", i }; state.editing = i; return; }
      }
    }
    drag = { mode: state.view3d ? "orbit" : "pan", sx, sy, yaw: state.yaw, tilt: state.tilt, C: state.camC.slice() };
    stage.classList.add("dragging");
  });
  stage.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const r = stage.getBoundingClientRect(), sx = e.clientX - r.left, sy = e.clientY - r.top;
    if (drag.mode === "orbit") {
      state.yaw = drag.yaw + (sx - drag.sx) * 0.008;
      state.tilt = clamp(drag.tilt + (sy - drag.sy) * 0.25, 0, 88);
    } else if (drag.mode === "pan") {
      const sc = scale();
      state.camC[0] = drag.C[0] - (sx - drag.sx) / sc; state.camC[1] = drag.C[1] + (sy - drag.sy) / sc;
      state.manualZoom = true; setOpt("autofit", false);
    } else if (drag.mode === "pos") {
      const w = unproject(sx, sy);
      state.disp[3 * drag.i] = w[0]; state.disp[3 * drag.i + 1] = w[1];
      drag.moved = true;
    } else if (drag.mode === "vel") {
      const p = project(state.disp[3 * drag.i], state.disp[3 * drag.i + 1], 0);
      state.dispV[3 * drag.i] = (sx - p[0]) / drag.vs;
      state.dispV[3 * drag.i + 1] = -((sy - p[1]) / drag.vs);
      drag.moved = true;
    }
  });
  const endDrag = () => {
    if (drag && (drag.mode === "pos" || drag.mode === "vel") && drag.moved) {
      setCustom(Array.from(state.disp), Array.from(state.dispV));
    }
    drag = null; state.editing = null; stage.classList.remove("dragging");
  };
  stage.addEventListener("pointerup", endDrag);
  stage.addEventListener("pointercancel", endDrag);
  stage.addEventListener("wheel", (e) => {
    e.preventDefault();
    state.camR *= Math.exp(e.deltaY * 0.0012);
    state.manualZoom = true; setOpt("autofit", false);
  }, { passive: false });
  stage.addEventListener("dblclick", () => { setOpt("autofit", true); });

  /* ================================================================ main loop */
  let last = performance.now(), frameN = 0, fpsAcc = 0, fpsN = 0;
  function frame(now) {
    const dt = Math.min(0.1, (now - last) / 1000);
    last = now; frameN++;
    if (state.playing && !state.collision) {
      const target = state.tDisp + state.speed * dt;
      const done = advanceTo(target, 14);
      state.lag = !done;
      state.tDisp = done ? target : state.sim.t;
    }
    if (drag && (drag.mode === "pos" || drag.mode === "vel")) { /* keep edited display state */ }
    else sampleDisplay();
    if (frameN % 2 === 0) pushDiag();
    if (frameN % 6 === 0 || !state.playing) classifyNow();
    updateCamera(dt);
    render();
    if (frameN % 3 === 0) {
      chartE.draw(state.diag.E); chartR.draw(state.diag.R); chartS.draw(state.twin ? state.diag.S : null); drawSphere();
      updateHud();
      fpsAcc += dt; fpsN++;
      $("perf").textContent = state.lag ? "solver-limited: slowed to keep accuracy" : "";
    }
    requestAnimationFrame(frame);
  }

  /* ================================================================ fate map */
  function setupFateMap() {
    if (!FATE) { $("fatemap").hidden = true; return; }
    const cv = $("fateCanvas"), ov = $("fateOverlay");
    const Wm = FATE.width, Hm = FATE.height;
    cv.width = Wm; cv.height = Hm;
    const bin = atob(FATE.data), data = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) data[i] = bin.charCodeAt(i);
    const none = [17, 22, 42];
    const total = data.length, cnt = FATE.counts;
    const pct = (n) => `${((100 * n) / total).toFixed(0)}%`;
    $("fateCount").textContent = total.toLocaleString();
    const RAMP = [[236, 240, 253], [120, 138, 233], [22, 27, 58]]; // fast → slow, one hue
    const paint = (mode) => {
      const img = cv.getContext("2d").createImageData(Wm, Hm);
      for (let i = 0; i < data.length; i++) {
        const b = data[i];
        let col;
        if (b === 255) col = mode === "when" ? [8, 10, 20] : none;
        else if (mode === "when") {
          const q = (b & 63) / 63, f = q < 0.5 ? q / 0.5 : (q - 0.5) / 0.5, A = q < 0.5 ? RAMP[0] : RAMP[1], B = q < 0.5 ? RAMP[1] : RAMP[2];
          col = [A[0] + (B[0] - A[0]) * f, A[1] + (B[1] - A[1]) * f, A[2] + (B[2] - A[2]) * f];
        } else {
          const k = b >> 6, q = (b & 63) / 63;
          const f = 1 - 0.8 * Math.pow(q, 0.8);
          const g = GLOW_RGB[k];
          col = [none[0] + (g[0] - none[0]) * f, none[1] + (g[1] - none[1]) * f, none[2] + (g[2] - none[2]) * f];
        }
        img.data[4 * i] = col[0]; img.data[4 * i + 1] = col[1]; img.data[4 * i + 2] = col[2]; img.data[4 * i + 3] = 255;
      }
      cv.getContext("2d").putImageData(img, 0, 0);
      $("fateLegend").innerHTML = mode === "when"
        ? `<span>escape soon<span class="ramp"></span>late (t = ${FATE.tLo} → ${FATE.tMax}, log scale)</span><span><span class="sw" style="background:rgb(8,10,20);box-shadow:inset 0 0 0 1px #2a3657"></span>no escape by t = ${FATE.tMax} (${pct(cnt.none)})</span>`
        : [0, 1, 2].map((i) => `<span><span class="sw" style="background:${GLOW[i]}"></span>body ${i + 1} escapes (${pct(cnt["body" + (i + 1)])})</span>`).join("") +
          `<span>bright = soon, dim = late</span><span><span class="sw" style="background:rgb(${none})"></span>no escape by t = ${FATE.tMax} (${pct(cnt.none)})</span>`;
      $("fateWhen").setAttribute("aria-pressed", String(mode === "when"));
      $("fateWho").setAttribute("aria-pressed", String(mode === "who"));
    };
    paint("when");
    $("fateWhen").addEventListener("click", () => paint("when"));
    $("fateWho").addEventListener("click", () => paint("who"));
    const toWorld = (c, r) => [FATE.xMin + ((FATE.xMax - FATE.xMin) * (c + 0.5)) / Wm, FATE.yMax * (1 - (r + 0.5) / Hm)];
    const drawOverlay = () => {
      const rr = cv.getBoundingClientRect(), w = rr.width, h = rr.height;
      ov.width = Math.round(w * DPR); ov.height = Math.round(h * DPR);
      const c = ov.getContext("2d"); c.setTransform(DPR, 0, 0, DPR, 0, 0); c.clearRect(0, 0, w, h);
      const X = (x) => ((x - FATE.xMin) / (FATE.xMax - FATE.xMin)) * w, Y = (y) => (1 - (y - FATE.yMin) / (FATE.yMax - FATE.yMin)) * h;
      c.font = "500 11px 'IBM Plex Mono', ui-monospace, monospace";
      for (const [i, x] of [[0, -0.5], [1, 0.5]]) {
        c.fillStyle = GLOW[i]; c.strokeStyle = "#060912"; c.lineWidth = 2;
        c.beginPath(); c.arc(X(x), Y(0) - 1, 6, 0, 2 * Math.PI); c.fill(); c.stroke();
        c.fillStyle = "#e6e9f2"; c.fillText(`body ${i + 1}`, X(x) - 22, Y(0) - 12);
      }
      const ly = Math.sqrt(3) / 2;
      c.strokeStyle = "rgba(230,233,242,0.9)"; c.lineWidth = 1.4;
      c.beginPath(); c.arc(X(0), Y(ly), 7, 0, 2 * Math.PI); c.stroke();
      c.fillStyle = "#e6e9f2"; c.fillText("Lagrange △", X(0) + 11, Y(ly) + 4);
      if (state.preset && state.preset.id === "freefall") {
        c.strokeStyle = "#ffffff"; c.lineWidth = 2;
        c.beginPath(); c.arc(X(state.ff[0]), Y(state.ff[1]), 8, 0, 2 * Math.PI); c.stroke();
      }
    };
    drawOverlay();
    window.addEventListener("resize", drawOverlay);
    const pick = (e) => {
      const rr = cv.getBoundingClientRect();
      const c = clamp(Math.floor(((e.clientX - rr.left) / rr.width) * Wm), 0, Wm - 1), r = clamp(Math.floor(((e.clientY - rr.top) / rr.height) * Hm), 0, Hm - 1);
      return { c, r, b: data[r * Wm + c], w: toWorld(c, r) };
    };
    cv.addEventListener("pointermove", (e) => {
      const p = pick(e);
      let fate = "no escape by t = " + FATE.tMax;
      if (p.b !== 255) { const k = p.b >> 6, q = p.b & 63; const t = FATE.tLo * Math.exp((q / 63) * Math.log(FATE.tMax / FATE.tLo)); fate = `<span style="color:${GLOW[k]}">●</span> body ${k + 1} escapes at t ≈ ${t.toFixed(t < 10 ? 1 : 0)}`; }
      showTip(e, `<div class="k">body 3 starts at (${p.w[0].toFixed(3)}, ${p.w[1].toFixed(3)})</div><div>${fate}</div><div class="k">click to run it</div>`);
    });
    cv.addEventListener("pointerleave", hideTip);
    cv.addEventListener("click", (e) => {
      const p = pick(e);
      state.ff = [p.w[0], p.w[1]];
      setOpt("twin", false);
      setOpt("autofit", true);
      state.eps = 0; $("eps").value = "0";
      loadPreset("freefall");
      state.playing = true; syncPlay();
      drawOverlay();
      hideTip();
      stage.scrollIntoView({ behavior: reduceMotion ? "auto" : "smooth", block: "center" });
    });
  }

  /* ================================================================ units converter */
  function setupConverter() {
    const G = 6.6743e-11, AU = 1.495978707e11, MSUN = 1.98847e30;
    const sys = $("cvSystem"), L = $("cvL"), M = $("cvM"), out = $("cvOut");
    const calc = () => {
      const l = +L.value * AU, mm = +M.value * MSUN;
      if (!(l > 0 && mm > 0)) { out.textContent = "Enter a positive length and mass."; return; }
      const ts = Math.sqrt((l * l * l) / (G * mm));
      const d = ts / 86400, y = d / 365.25;
      const human = y >= 1 ? `${fmt(y, y < 100 ? 2 : 0)} years` : d >= 1 ? `${fmt(d, 2)} days` : `${fmt(ts / 3600, 2)} hours`;
      out.textContent = `t* = ${human}  ·  one simulator time unit. A circular orbit of radius L takes 2π t* = ${y * 2 * Math.PI >= 1 ? fmt(2 * Math.PI * y, 2) + " years" : fmt(2 * Math.PI * d, 2) + " days"}.`;
    };
    sys.addEventListener("change", () => { const [, l, m] = sys.value.split(",").map(Number); L.value = String(l); M.value = String(m); calc(); });
    L.addEventListener("input", calc); M.addEventListener("input", calc);
    $("converter").addEventListener("submit", (e) => e.preventDefault());
    calc();
  }

  /* ================================================================ boot */
  function boot() {
    resize();
    window.addEventListener("resize", resize);
    if (window.ResizeObserver) new ResizeObserver(resize).observe(stage);
    populatePresetSelect();
    loadPreset("figure8");
    syncPlay();
    setupFateMap();
    setupConverter();
    requestAnimationFrame((t) => { last = t; frame(t); });
  }
  window.__threeBody = { state, PRESETS, loadPreset, setOpt, restart, dims: () => ({ W, H, DPR, scale: scale(), p3: project(state.disp[6], state.disp[7], 0) }) };
  boot();
})();
