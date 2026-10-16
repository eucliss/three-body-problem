/* The Three-Body Ladder: findings page. Every orbit is drawn from its own integrated path
 * (build_data.js), the words are spelled at the syzygy times measured on that path, and the
 * Check buttons re-integrate an orbit live with the Observatory's solver. */
(function () {
  "use strict";
  const D = window.LADDER;
  const $ = (id) => document.getElementById(id);
  const REDUCED = matchMedia("(prefers-reduced-motion: reduce)").matches;
  const DPR = () => Math.min(window.devicePixelRatio || 1, 2);
  const C = {
    body: ["#c2841a", "#3a91d4", "#de5577"], glow: ["#f6c66a", "#8fd0ff", "#ff9bb3"],
    rgb: [[194, 132, 26], [58, 145, 212], [222, 85, 119]], grgb: [[246, 198, 106], [143, 208, 255], [255, 155, 179]],
    neu: "#8ef0cf", known: "#8590ad", ink: "#ecebe6", soft: "#cdd2df", muted: "#959db4", faint: "#5c6582",
    line: "#172038", lineStrong: "#26314f", accent: "#aab8ff", pink: "#ff9bb3",
  };
  const RATE = 2.6; // orbit time units per second at 1x
  const TAIL = 4.6; // comet tail length, in orbit time units
  const rgba = (c, a) => `rgba(${c[0]},${c[1]},${c[2]},${a})`;
  const fmt = (x, n) => (+x).toFixed(n);
  const mod = (a, b) => ((a % b) + b) % b;

  /* ---------------- data ---------------- */
  function b64(s, Type) {
    const bin = atob(s), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return new Type(u8.buffer);
  }
  function shapeOf(p, i) {
    const o = i * 6, x1 = p[o], y1 = p[o + 1], x2 = p[o + 2], y2 = p[o + 3], x3 = p[o + 4], y3 = p[o + 5];
    const r0 = (x1 - x2) / Math.SQRT2, r1 = (y1 - y2) / Math.SQRT2, s6 = Math.sqrt(6);
    const l0 = (x1 + x2 - 2 * x3) / s6, l1 = (y1 + y2 - 2 * y3) / s6;
    const R2 = r0 * r0 + r1 * r1 + l0 * l0 + l1 * l1;
    return [(2 * (r0 * l0 + r1 * l1)) / R2, (l0 * l0 + l1 * l1 - r0 * r0 - r1 * r1) / R2, (2 * (r0 * l1 - r1 * l0)) / R2];
  }
  const ORB = {};
  for (const o of D.orbits) {
    o.ts = b64(o.t, Float32Array);
    const q = b64(o.xy, Int16Array), [cx, cy, span] = o.box;
    o.P = new Float32Array(q.length);
    for (let i = 0; i < q.length; i++) o.P[i] = (q[i] / 65520) * span + (i & 1 ? cy : cx);
    o.N = o.ts.length;
    o.S = new Float32Array(o.N * 3);
    for (let i = 0; i < o.N; i++) o.S.set(shapeOf(o.P, i), i * 3);
    o.t0 = 0; // loops start at the starting line-up, where the pair has not turned yet
    o.tauL = o.spell.map(([t]) => t);
    // turns of the heavy pair: unwrapped angle of r2 - r1 over 2 pi
    o.TH = new Float32Array(o.N);
    let acc = 0, prev = Math.atan2(o.P[3] - o.P[1], o.P[2] - o.P[0]);
    for (let i = 1; i < o.N; i++) {
      const a = Math.atan2(o.P[i * 6 + 3] - o.P[i * 6 + 1], o.P[i * 6 + 2] - o.P[i * 6]);
      let d = a - prev; while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI;
      acc += d; prev = a; o.TH[i] = acc / (2 * Math.PI);
    }
    o.turnsMax = Math.max(...o.TH.map(Math.abs));
    o.label = o.extra ? "3′" : String(o.k);
    ORB[o.id] = o;
  }
  const MASSES = [...new Set(D.orbits.map((o) => o.m3))].filter((m) => D.orbits.filter((o) => o.m3 === m).length > 1).sort((a, b) => a - b);
  const byMass = (m3) => D.orbits.filter((o) => Math.abs(o.m3 - m3) < 1e-9).sort((a, b) => a.k - b.k || (a.extra ? 1 : -1));

  /* path position at orbit time t (linear between stored samples); returns the sample index */
  function posAt(o, t, out) {
    const ts = o.ts, N = o.N;
    t = mod(t, o.T);
    let lo = 0, hi = N - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (ts[m] <= t) lo = m; else hi = m; }
    const f = (t - ts[lo]) / (ts[hi] - ts[lo] || 1);
    for (let c = 0; c < 6; c++) out[c] = o.P[lo * 6 + c] + f * (o.P[hi * 6 + c] - o.P[lo * 6 + c]);
    return lo;
  }
  function turnsAt(o, t) {
    const ts = o.ts; t = mod(t, o.T);
    let lo = 0, hi = o.N - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (ts[m] <= t) lo = m; else hi = m; }
    const f = (t - ts[lo]) / (ts[hi] - ts[lo] || 1);
    return o.TH[lo] + f * (o.TH[hi] - o.TH[lo]);
  }
  /* times and positions of the path over [t - span, t], following the period around */
  function window6(o, t, span) {
    const T = o.T, ts = o.ts, N = o.N, outT = [], outP = [];
    const a = t - span, p = new Float64Array(6);
    posAt(o, a, p); outT.push(a); outP.push(...p);
    let base = Math.floor(a / T) * T, j;
    { const tm = a - base; let lo = 0, hi = N - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (ts[m] <= tm) lo = m; else hi = m; } j = hi; }
    let guard = 0;
    while (guard++ < 200000) {
      if (j >= N - 1) { j = 1; base += T; }
      const tj = base + ts[j];
      if (tj >= t) break;
      outT.push(tj); for (let c = 0; c < 6; c++) outP.push(o.P[j * 6 + c]);
      j++;
    }
    posAt(o, t, p); outT.push(t); outP.push(...p);
    return { T: outT, P: outP };
  }

  /* ---------------- canvas helpers ---------------- */
  function fit(cv) {
    const r = cv.getBoundingClientRect(), d = DPR();
    const w = Math.max(1, Math.round(r.width * d)), h = Math.max(1, Math.round(r.height * d));
    if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; return true; }
    return false;
  }
  function viewFor(o, w, h, fill, dx, dy) {
    const s = (Math.min(w, h) * fill) / o.box[2];
    return { s, X: (x) => w / 2 + (dx || 0) + (x - o.box[0]) * s, Y: (y) => h / 2 + (dy || 0) - (y - o.box[1]) * s };
  }
  function ghost(ctx, o, v, d, alpha, wide) {
    ctx.save(); ctx.globalCompositeOperation = "lighter"; ctx.lineJoin = ctx.lineCap = "round";
    for (let b = 0; b < 3; b++) {
      const passes = wide ? [[7, 0.03], [2.4, 0.06], [1, alpha]] : [[2.2, 0.08 * alpha * 3], [0.9, alpha]];
      for (const [lw, a] of passes) {
        ctx.beginPath();
        for (let i = 0; i < o.N; i++) {
          const X = v.X(o.P[i * 6 + 2 * b]), Y = v.Y(o.P[i * 6 + 2 * b + 1]);
          i ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y);
        }
        ctx.strokeStyle = rgba(b === 2 ? C.grgb[2] : C.grgb[b], a); ctx.lineWidth = lw * d; ctx.stroke();
      }
    }
    ctx.restore();
  }
  function tails(ctx, o, v, t, span, d, wScale) {
    const W = window6(o, t, span), n = W.T.length;
    if (n < 2) return;
    ctx.save(); ctx.globalCompositeOperation = "lighter"; ctx.lineJoin = ctx.lineCap = "round";
    const NB = 14;
    for (let b = 0; b < 3; b++) {
      for (let pass = 0; pass < 2; pass++) {
        for (let k = 0; k < NB; k++) {
          const f0 = k / NB, f1 = (k + 1) / NB;
          ctx.beginPath(); let started = false;
          for (let i = 0; i < n; i++) {
            const age = (W.T[n - 1] - W.T[i]) / span, fr = 1 - age;
            if (fr < f0 - 0.0001 && (i + 1 >= n || 1 - (W.T[n - 1] - W.T[i + 1]) / span < f0)) continue;
            if (fr > f1 + 1e-9 && started) { const X = v.X(W.P[i * 6 + 2 * b]), Y = v.Y(W.P[i * 6 + 2 * b + 1]); ctx.lineTo(X, Y); break; }
            const X = v.X(W.P[i * 6 + 2 * b]), Y = v.Y(W.P[i * 6 + 2 * b + 1]);
            if (!started) { ctx.moveTo(X, Y); started = true; } else ctx.lineTo(X, Y);
          }
          if (!started) continue;
          const fm = (f0 + f1) / 2, ease = fm * fm;
          if (pass === 0) { ctx.strokeStyle = rgba(C.grgb[b], 0.11 * ease); ctx.lineWidth = (2 + 5 * fm) * d * wScale; }
          else { ctx.strokeStyle = rgba(fm > 0.85 ? [255, 250, 240] : C.grgb[b], 0.12 + 0.88 * ease); ctx.lineWidth = (0.5 + 1.5 * fm) * d * wScale; }
          ctx.stroke();
        }
      }
    }
    ctx.restore();
  }
  function bodies(ctx, v, p, m3, d, size) {
    ctx.save(); ctx.globalCompositeOperation = "lighter";
    for (let b = 0; b < 3; b++) {
      const X = v.X(p[2 * b]), Y = v.Y(p[2 * b + 1]);
      const r = size * d * (b < 2 ? 1 : Math.max(0.42, Math.cbrt(m3)));
      const g = ctx.createRadialGradient(X, Y, r * 0.6, X, Y, r * 4.2);
      g.addColorStop(0, rgba(C.grgb[b], 0.32)); g.addColorStop(0.35, rgba(C.grgb[b], 0.08)); g.addColorStop(1, rgba(C.grgb[b], 0));
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(X, Y, r * 4.2, 0, 7); ctx.fill();
    }
    ctx.restore();
    for (let b = 0; b < 3; b++) {
      const X = v.X(p[2 * b]), Y = v.Y(p[2 * b + 1]);
      const r = size * d * (b < 2 ? 1 : Math.max(0.42, Math.cbrt(m3)));
      const g = ctx.createRadialGradient(X - r * 0.3, Y - r * 0.3, 0, X, Y, r);
      g.addColorStop(0, "#ffffff"); g.addColorStop(0.45, C.glow[b]); g.addColorStop(1, C.body[b]);
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(X, Y, r, 0, 7); ctx.fill();
    }
  }
  /* line-ups: a fading line through the three stars at each syzygy, and the letter it writes */
  function syzygyFx(ctx, o, v, t, d, withLetters) {
    const dur = 1.5, p = new Float64Array(6);
    ctx.save(); ctx.globalCompositeOperation = "lighter"; ctx.lineCap = "round";
    for (const [ts, body] of o.syzygies) {
      const age = mod(t - ts, o.T);
      if (age > dur) continue;
      const a = 1 - age / dur;
      posAt(o, ts, p);
      const xs = [0, 1, 2].map((b) => [v.X(p[2 * b]), v.Y(p[2 * b + 1])]);
      let i0 = 0, i1 = 1, best = -1;
      for (const [i, j] of [[0, 1], [0, 2], [1, 2]]) { const q = Math.hypot(xs[i][0] - xs[j][0], xs[i][1] - xs[j][1]); if (q > best) { best = q; i0 = i; i1 = j; } }
      const ex = (xs[i1][0] - xs[i0][0]) * 0.22, ey = (xs[i1][1] - xs[i0][1]) * 0.22;
      const g = ctx.createLinearGradient(xs[i0][0] - ex, xs[i0][1] - ey, xs[i1][0] + ex, xs[i1][1] + ey);
      const c = C.grgb[body - 1];
      g.addColorStop(0, rgba(c, 0)); g.addColorStop(0.2, rgba(c, 0.55 * a)); g.addColorStop(0.5, rgba([255, 255, 255], 0.5 * a)); g.addColorStop(0.8, rgba(c, 0.55 * a)); g.addColorStop(1, rgba(c, 0));
      ctx.strokeStyle = g; ctx.lineWidth = 1.2 * d;
      ctx.beginPath(); ctx.moveTo(xs[i0][0] - ex, xs[i0][1] - ey); ctx.lineTo(xs[i1][0] + ex, xs[i1][1] + ey); ctx.stroke();
      const m = xs[body - 1], rr = (6 + 26 * (1 - a)) * d;
      ctx.strokeStyle = rgba(c, 0.6 * a); ctx.lineWidth = 1 * d; ctx.beginPath(); ctx.arc(m[0], m[1], rr, 0, 7); ctx.stroke();
    }
    ctx.restore();
    if (!withLetters) return;
    ctx.save(); ctx.textAlign = "center"; ctx.textBaseline = "middle";
    const seen = new Map();
    for (const [ts, ch] of o.spell) {
      const age = mod(t - ts, o.T);
      if (age > 2.2) continue;
      const a = Math.min(1, (2.2 - age) / 1.2), n = seen.get(ts) || 0; seen.set(ts, n + 1);
      posAt(o, ts, p);
      const X = v.X((p[0] + p[2] + p[4]) / 3) + n * 15 * d, Y = v.Y((p[1] + p[3] + p[5]) / 3) - (22 + 30 * (age / 2.2)) * d;
      ctx.font = `600 ${Math.round(19 * d)}px "IBM Plex Mono", ui-monospace, monospace`;
      ctx.shadowColor = "rgba(143,208,255,0.9)"; ctx.shadowBlur = 14 * d;
      ctx.fillStyle = `rgba(255,255,255,${0.95 * a})`; ctx.fillText(ch, X, Y);
    }
    ctx.restore();
  }

  /* ---------------- sky ---------------- */
  function rng(seed) { let s = seed >>> 0; return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }; }
  function drawSky(cv) {
    fit(cv);
    const ctx = cv.getContext("2d"), w = cv.width, h = cv.height, d = DPR(), R = rng(7);
    ctx.clearRect(0, 0, w, h);
    const neb = [[0.74, 0.28, 0.62, [92, 78, 190], 0.15], [0.18, 0.78, 0.55, [40, 140, 150], 0.09], [0.52, 0.62, 0.35, [200, 70, 120], 0.05], [0.3, 0.2, 0.4, [60, 90, 200], 0.06]];
    for (const [x, y, r, c, a] of neb) {
      const g = ctx.createRadialGradient(w * x, h * y, 0, w * x, h * y, Math.max(w, h) * r);
      g.addColorStop(0, rgba(c, a)); g.addColorStop(0.5, rgba(c, a * 0.35)); g.addColorStop(1, rgba(c, 0));
      ctx.fillStyle = g; ctx.fillRect(0, 0, w, h);
    }
    const n = Math.round((w * h) / (2100 * d * d));
    for (let i = 0; i < n; i++) {
      let x = R() * w, y = R() * h;
      if (R() < 0.45) { const u = R(); x = u * w; y = h * (0.92 - 0.8 * u) + (R() + R() + R() - 1.5) * h * 0.09; } // a faint band
      const mag = Math.pow(R(), 3.2), r = (0.25 + 1.25 * mag) * d, temp = R();
      const col = temp < 0.2 ? [255, 214, 170] : temp < 0.4 ? [190, 210, 255] : [235, 238, 255];
      ctx.fillStyle = rgba(col, 0.18 + 0.75 * mag);
      ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.fill();
      if (mag > 0.82) {
        const g = ctx.createRadialGradient(x, y, 0, x, y, r * 6);
        g.addColorStop(0, rgba(col, 0.25)); g.addColorStop(1, rgba(col, 0));
        ctx.fillStyle = g; ctx.beginPath(); ctx.arc(x, y, r * 6, 0, 7); ctx.fill();
      }
    }
    const vg = ctx.createRadialGradient(w / 2, h / 2, Math.min(w, h) * 0.35, w / 2, h / 2, Math.max(w, h) * 0.75);
    vg.addColorStop(0, "rgba(3,5,12,0)"); vg.addColorStop(1, "rgba(3,5,12,0.75)");
    ctx.fillStyle = vg; ctx.fillRect(0, 0, w, h);
  }

  /* ---------------- shape sphere ---------------- */
  const PUNC = [[[0, 1, 0], "1·2"], [[-Math.sqrt(3) / 2, -0.5, 0], "1·3"], [[Math.sqrt(3) / 2, -0.5, 0], "2·3"]];
  function rotator(yaw, pitch) {
    const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
    return (n0, n1, n2) => { // world: right = n0, up = n2 (the equilateral poles), depth = n1
      const x1 = cy * n0 + sy * n1, z1 = -sy * n0 + cy * n1;
      return [x1, cp * n2 - sp * z1, sp * n2 + cp * z1];
    };
  }
  function drawSphere(cv, o, t, yaw, pitch, opts) {
    const ctx = cv.getContext("2d"), w = cv.width, h = cv.height, d = DPR();
    ctx.clearRect(0, 0, w, h);
    const R = Math.min(w, h) * (opts.r || 0.4), cx = w / 2, cy = h / 2 + (opts.dy || 0) * h, rot = rotator(yaw, pitch);
    const P = (n0, n1, n2) => { const q = rot(n0, n1, n2); return [cx + q[0] * R, cy - q[1] * R, q[2]]; };
    // body
    const g = ctx.createRadialGradient(cx - R * 0.35, cy - R * 0.42, R * 0.05, cx, cy, R);
    g.addColorStop(0, "rgba(150,165,255,0.20)"); g.addColorStop(0.55, "rgba(40,55,120,0.20)"); g.addColorStop(1, "rgba(8,12,30,0.55)");
    ctx.fillStyle = g; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.fill();
    const halo = ctx.createRadialGradient(cx, cy, R * 0.96, cx, cy, R * 1.25);
    halo.addColorStop(0, "rgba(170,184,255,0.16)"); halo.addColorStop(1, "rgba(170,184,255,0)");
    ctx.fillStyle = halo; ctx.beginPath(); ctx.arc(cx, cy, R * 1.25, 0, 7); ctx.fill();
    ctx.strokeStyle = "rgba(170,184,255,0.28)"; ctx.lineWidth = 1 * d; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.stroke();
    // graticule
    const seg = (pts, front, back, lw) => {
      for (let i = 1; i < pts.length; i++) {
        const a = pts[i - 1], b = pts[i], z = (a[2] + b[2]) / 2;
        ctx.strokeStyle = z >= 0 ? front : back; ctx.lineWidth = lw * d;
        ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
      }
    };
    for (const lat of [-60, -30, 30, 60]) {
      const pts = [], c = Math.cos((lat * Math.PI) / 180), s = Math.sin((lat * Math.PI) / 180);
      for (let i = 0; i <= 96; i++) { const f = (i / 96) * 2 * Math.PI; pts.push(P(c * Math.cos(f), c * Math.sin(f), s)); }
      seg(pts, "rgba(170,184,255,0.13)", "rgba(170,184,255,0.04)", 0.8);
    }
    for (let k = 0; k < 6; k++) {
      const pts = [], f = (k / 6) * Math.PI;
      for (let i = 0; i <= 96; i++) { const th = (i / 96) * 2 * Math.PI; pts.push(P(Math.cos(th) * Math.cos(f), Math.cos(th) * Math.sin(f), Math.sin(th))); }
      seg(pts, "rgba(170,184,255,0.10)", "rgba(170,184,255,0.035)", 0.8);
    }
    // equator: three arcs, coloured by the star in the middle of the line-up
    const arcs = [[90, 210, 0], [210, 330, 2], [330, 450, 1]];
    for (const [a0, a1, b] of arcs) {
      const pts = [];
      for (let i = 0; i <= 60; i++) { const f = ((a0 + ((a1 - a0) * i) / 60) * Math.PI) / 180; pts.push(P(Math.cos(f), Math.sin(f), 0)); }
      seg(pts, rgba(C.grgb[b], 0.85), rgba(C.grgb[b], 0.22), opts.big ? 2.2 : 1.6);
    }
    // the path on the sphere: back half dim, front half bright, recent part as a comet
    ctx.lineJoin = ctx.lineCap = "round";
    const S = o.S, N = o.N;
    let prev = P(S[0], S[1], S[2]);
    for (let i = 1; i < N; i++) {
      const q = P(S[i * 3], S[i * 3 + 1], S[i * 3 + 2]), z = (q[2] + prev[2]) / 2;
      ctx.strokeStyle = z >= 0 ? "rgba(225,232,255,0.42)" : "rgba(225,232,255,0.09)"; ctx.lineWidth = (z >= 0 ? 1.1 : 0.8) * d;
      ctx.beginPath(); ctx.moveTo(prev[0], prev[1]); ctx.lineTo(q[0], q[1]); ctx.stroke();
      prev = q;
    }
    const W = window6(o, t, TAIL * 0.9), n = W.T.length;
    const sp = [];
    for (let i = 0; i < n; i++) { const s = shapeOf(W.P, i); sp.push(P(s[0], s[1], s[2])); }
    ctx.save(); ctx.globalCompositeOperation = "lighter";
    for (let i = 1; i < n; i++) {
      const fr = 1 - (W.T[n - 1] - W.T[i]) / (TAIL * 0.9), z = sp[i][2];
      ctx.strokeStyle = `rgba(160,220,255,${(z >= 0 ? 0.9 : 0.35) * fr * fr})`; ctx.lineWidth = (0.8 + 2.4 * fr) * d;
      ctx.beginPath(); ctx.moveTo(sp[i - 1][0], sp[i - 1][1]); ctx.lineTo(sp[i][0], sp[i][1]); ctx.stroke();
    }
    // crossing flashes
    for (const [ts, body] of o.syzygies) {
      const age = mod(t - ts, o.T);
      if (age > 1.4) continue;
      const a = 1 - age / 1.4, p6 = new Float64Array(6); posAt(o, ts, p6);
      const s = shapeOf(p6, 0), q = P(s[0], s[1], s[2]);
      const gg = ctx.createRadialGradient(q[0], q[1], 0, q[0], q[1], (6 + 18 * (1 - a)) * d);
      gg.addColorStop(0, rgba(C.grgb[body - 1], 0.9 * a)); gg.addColorStop(1, rgba(C.grgb[body - 1], 0));
      ctx.fillStyle = gg; ctx.beginPath(); ctx.arc(q[0], q[1], (6 + 18 * (1 - a)) * d, 0, 7); ctx.fill();
    }
    const head = sp[n - 1];
    const hg = ctx.createRadialGradient(head[0], head[1], 0, head[0], head[1], 9 * d);
    hg.addColorStop(0, "rgba(255,255,255,0.95)"); hg.addColorStop(0.3, "rgba(160,220,255,0.6)"); hg.addColorStop(1, "rgba(160,220,255,0)");
    ctx.fillStyle = hg; ctx.beginPath(); ctx.arc(head[0], head[1], 9 * d, 0, 7); ctx.fill();
    ctx.restore();
    // collision points
    ctx.font = `500 ${Math.round((opts.big ? 11 : 9) * d)}px "IBM Plex Mono", monospace`; ctx.textAlign = "center"; ctx.textBaseline = "middle";
    for (const [n3, lab] of PUNC) {
      const q = P(n3[0], n3[1], n3[2]), front = q[2] >= 0;
      ctx.fillStyle = "#02030a"; ctx.strokeStyle = front ? "rgba(236,235,230,0.85)" : "rgba(236,235,230,0.25)"; ctx.lineWidth = 1.3 * d;
      ctx.beginPath(); ctx.arc(q[0], q[1], (opts.big ? 4.5 : 3.2) * d, 0, 7); ctx.fill(); ctx.stroke();
      if (opts.big && front) { ctx.fillStyle = "rgba(205,210,223,0.8)"; ctx.fillText(lab, q[0], q[1] - 13 * d); }
    }
    if (opts.big) {
      for (const [z, lab] of [[1, "equilateral"], [-1, "equilateral, mirrored"]]) {
        const q = P(0, 0, z);
        if (q[2] < -0.2) continue;
        ctx.fillStyle = "rgba(170,184,255,0.55)"; ctx.beginPath(); ctx.arc(q[0], q[1], 2 * d, 0, 7); ctx.fill();
        ctx.fillStyle = "rgba(149,157,180,0.75)"; ctx.fillText(lab, q[0], q[1] + (z > 0 ? -12 : 12) * d);
      }
    }
  }

  /* ---------------- hero ---------------- */
  const stage = $("stage"), sky = $("sky"), ghostCv = $("ghost"), live = $("live"), lens = $("lens"), spinCv = $("spin");
  let mass = 0.5, cur = ORB["L4@0.5"] || D.orbits[0], tau = 0, playing = !REDUCED, speed = 1;
  let dragging = false, fYaw = 0.9, fPitch = 0.38;
  const isNarrow = () => matchMedia("(max-width: 1100px)").matches;

  function heroView() {
    fit(live);
    const w = live.width, h = live.height, narrow = w < 900 * DPR();
    return viewFor(cur, w, h, narrow ? 0.8 : 0.74, narrow ? 0 : -w * 0.03, narrow ? h * 0.02 : h * 0.01);
  }
  function renderGhost() {
    fit(ghostCv);
    const ctx = ghostCv.getContext("2d");
    ctx.clearRect(0, 0, ghostCv.width, ghostCv.height);
    ghost(ctx, cur, heroView(), DPR(), 0.3, true);
  }
  function renderHero() {
    const ctx = live.getContext("2d"), d = DPR();
    ctx.clearRect(0, 0, live.width, live.height);
    const v = heroView(), t = cur.t0 + tau, p = new Float64Array(6);
    syzygyFx(ctx, cur, v, t, d, true);
    tails(ctx, cur, v, t, TAIL, d, 1);
    posAt(cur, t, p);
    bodies(ctx, v, p, cur.m3, d, Math.max(3.4, Math.min(6, v.s * 0.075 / d)));
    $("clock").innerHTML = `t <b>${fmt(tau, 2)}</b> of <b>${fmt(cur.T, 2)}</b>`;
    if (lens.offsetParent) drawLens(t, p);
    drawSpin(t);
    spellUpdate();
  }
  /* the lens: the heavy pair in its own centre-of-mass frame, with a turn counter */
  let lensCache = null, lensTxt = "";
  function drawLens(t, p) {
    fit(lens);
    const ctx = lens.getContext("2d"), w = lens.width, h = lens.height, d = DPR(), R = w * 0.36, cx = w / 2, cy = h / 2;
    if (!lensCache || lensCache.o !== cur || lensCache.w !== w) {
      lensCache = { o: cur, w, cv: document.createElement("canvas") };
      lensCache.cv.width = w; lensCache.cv.height = h;
      const g = lensCache.cv.getContext("2d");
      g.globalCompositeOperation = "lighter"; g.lineJoin = g.lineCap = "round";
      for (const b of [0, 1]) {
        g.beginPath();
        for (let i = 0; i < cur.N; i++) {
          const mx = (cur.P[i * 6] + cur.P[i * 6 + 2]) / 2, my = (cur.P[i * 6 + 1] + cur.P[i * 6 + 3]) / 2;
          const X = cx + (cur.P[i * 6 + 2 * b] - mx) * R, Y = cy - (cur.P[i * 6 + 2 * b + 1] - my) * R;
          i ? g.lineTo(X, Y) : g.moveTo(X, Y);
        }
        g.strokeStyle = rgba(C.grgb[b], 0.16); g.lineWidth = 0.9 * d; g.stroke();
      }
    }
    ctx.clearRect(0, 0, w, h);
    ctx.drawImage(lensCache.cv, 0, 0);
    const mx = (p[0] + p[2]) / 2, my = (p[1] + p[3]) / 2;
    // direction to the light star, on the rim
    const a3 = Math.atan2(p[5] - my, p[4] - mx);
    ctx.strokeStyle = rgba(C.grgb[2], 0.9); ctx.lineWidth = 2.2 * d; ctx.lineCap = "round";
    ctx.beginPath(); ctx.arc(cx, cy, w * 0.46, -a3 - 0.16, -a3 + 0.16); ctx.stroke();
    // turn gauge: an arc that winds once per turn
    const tn = turnsAt(cur, t), frac = tn - Math.trunc(tn);
    ctx.strokeStyle = "rgba(170,184,255,0.12)"; ctx.lineWidth = 1.5 * d;
    ctx.beginPath(); ctx.arc(cx, cy, w * 0.42, 0, 7); ctx.stroke();
    ctx.strokeStyle = "rgba(170,184,255,0.75)"; ctx.lineWidth = 2 * d;
    ctx.beginPath(); ctx.arc(cx, cy, w * 0.42, 0, -2 * Math.PI * frac, frac > 0); ctx.stroke();
    // the pair, with a faint bar between them
    const pts = [0, 1].map((b) => [cx + (p[2 * b] - mx) * R, cy - (p[2 * b + 1] - my) * R]);
    ctx.strokeStyle = "rgba(236,235,230,0.18)"; ctx.lineWidth = 1 * d;
    ctx.beginPath(); ctx.moveTo(pts[0][0], pts[0][1]); ctx.lineTo(pts[1][0], pts[1][1]); ctx.stroke();
    const v = { X: (x) => cx + (x - mx) * R, Y: (y) => cy - (y - my) * R, s: R };
    tailsLocal(ctx, t, mx, my, R, cx, cy, d);
    for (const b of [0, 1]) {
      const r = 4.2 * d, [X, Y] = pts[b];
      const g = ctx.createRadialGradient(X - r * 0.3, Y - r * 0.3, 0, X, Y, r);
      g.addColorStop(0, "#fff"); g.addColorStop(0.45, C.glow[b]); g.addColorStop(1, C.body[b]);
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(X, Y, r, 0, 7); ctx.fill();
    }
    const txt = `${Math.abs(tn).toFixed(1)} of ${cur.k} turn${cur.k === 1 ? "" : "s"}`;
    if (txt !== lensTxt) { lensTxt = txt; $("lensTurns").textContent = txt; }
    void v;
  }
  function tailsLocal(ctx, t, mx0, my0, R, cx, cy, d) { // comet tails of the pair in the lens frame
    const W = window6(cur, t, TAIL * 0.7), n = W.T.length;
    ctx.save(); ctx.globalCompositeOperation = "lighter"; ctx.lineCap = "round";
    for (const b of [0, 1]) for (let i = 1; i < n; i++) {
      const fr = 1 - (W.T[n - 1] - W.T[i]) / (TAIL * 0.7);
      const m0x = (W.P[(i - 1) * 6] + W.P[(i - 1) * 6 + 2]) / 2, m0y = (W.P[(i - 1) * 6 + 1] + W.P[(i - 1) * 6 + 3]) / 2;
      const m1x = (W.P[i * 6] + W.P[i * 6 + 2]) / 2, m1y = (W.P[i * 6 + 1] + W.P[i * 6 + 3]) / 2;
      ctx.strokeStyle = rgba(C.grgb[b], 0.85 * fr * fr); ctx.lineWidth = (0.6 + 1.6 * fr) * d;
      ctx.beginPath();
      ctx.moveTo(cx + (W.P[(i - 1) * 6 + 2 * b] - m0x) * R, cy - (W.P[(i - 1) * 6 + 2 * b + 1] - m0y) * R);
      ctx.lineTo(cx + (W.P[i * 6 + 2 * b] - m1x) * R, cy - (W.P[i * 6 + 2 * b + 1] - m1y) * R);
      ctx.stroke();
    }
    ctx.restore();
  }
  /* spin track: turns of the pair through the period, with the letters at their times */
  let spinCache = null;
  function drawSpin(t) {
    if (!spinCv.offsetParent) return;
    const rz = fit(spinCv), w = spinCv.width, h = spinCv.height, d = DPR();
    const padL = 4 * d, padR = 4 * d, top = 6 * d, bot = 18 * d;
    const X = (tt) => padL + (tt / cur.T) * (w - padL - padR), sgn = Math.sign(cur.TH[Math.floor(cur.N / 2)]) || 1;
    const Y = (q) => top + (1 - (sgn * q) / Math.max(1, cur.turnsMax)) * (h - top - bot);
    if (rz || !spinCache || spinCache.o !== cur) {
      spinCache = { o: cur, cv: document.createElement("canvas") };
      spinCache.cv.width = w; spinCache.cv.height = h;
      const g = spinCache.cv.getContext("2d");
      g.strokeStyle = "rgba(38,49,79,0.9)"; g.lineWidth = 1;
      for (let q = 0; q <= Math.round(cur.turnsMax); q++) { g.beginPath(); g.moveTo(padL, Y(sgn * q)); g.lineTo(w - padR, Y(sgn * q)); g.stroke(); }
      const grad = g.createLinearGradient(0, top, 0, h - bot);
      grad.addColorStop(0, "rgba(170,184,255,0.30)"); grad.addColorStop(1, "rgba(170,184,255,0.02)");
      g.beginPath(); g.moveTo(X(0), Y(0));
      for (let i = 0; i < cur.N; i++) g.lineTo(X(cur.ts[i]), Y(cur.TH[i]));
      g.lineTo(X(cur.T), Y(0)); g.closePath(); g.fillStyle = grad; g.fill();
      g.beginPath();
      for (let i = 0; i < cur.N; i++) i ? g.lineTo(X(cur.ts[i]), Y(cur.TH[i])) : g.moveTo(X(cur.ts[i]), Y(cur.TH[i]));
      g.strokeStyle = "rgba(170,184,255,0.85)"; g.lineWidth = 1.4 * d; g.stroke();
      g.font = `500 ${Math.round(10 * d)}px "IBM Plex Mono", monospace`; g.textAlign = "center"; g.textBaseline = "alphabetic";
      const seen = new Map();
      for (const [ts, ch] of cur.spell) {
        const n = seen.get(ts) || 0; seen.set(ts, n + 1);
        g.fillStyle = ch === ch.toUpperCase() ? "rgba(205,210,223,0.75)" : "rgba(149,157,180,0.75)";
        g.fillText(ch, Math.min(w - 5 * d, Math.max(5 * d, X(ts) + (n ? 7 * d : 0))), h - 4 * d);
      }
      g.fillStyle = "rgba(149,157,180,0.7)"; g.textAlign = "left";
      g.fillText(`${Math.round(cur.turnsMax)} turn${Math.round(cur.turnsMax) === 1 ? "" : "s"}`, padL + 2 * d, top + 9 * d);
      g.fillText("0", padL + 2 * d, Y(0) - 3 * d);
      g.strokeStyle = "rgba(149,157,180,0.25)"; g.setLineDash([2 * d, 3 * d]); g.beginPath(); g.moveTo(X(cur.T / 2), top); g.lineTo(X(cur.T / 2), h - bot); g.stroke(); g.setLineDash([]);
    }
    const ctx = spinCv.getContext("2d");
    ctx.clearRect(0, 0, w, h); ctx.drawImage(spinCache.cv, 0, 0);
    const tt = mod(t, cur.T), x = X(tt), y = Y(turnsAt(cur, tt));
    ctx.strokeStyle = "rgba(236,235,230,0.35)"; ctx.lineWidth = 1 * d;
    ctx.beginPath(); ctx.moveTo(x, top - 2 * d); ctx.lineTo(x, h - bot + 2 * d); ctx.stroke();
    const g2 = ctx.createRadialGradient(x, y, 0, x, y, 8 * d);
    g2.addColorStop(0, "rgba(255,255,255,0.95)"); g2.addColorStop(0.35, "rgba(170,184,255,0.7)"); g2.addColorStop(1, "rgba(170,184,255,0)");
    ctx.fillStyle = g2; ctx.beginPath(); ctx.arc(x, y, 8 * d, 0, 7); ctx.fill();
  }

  /* word bar */
  let wordSpans = [], lastLit = -1;
  function groupsOf(w) {
    const m = /^((?:Ba)+)((?:bA)+)$/.exec(w);
    if (!m || m[1].length !== m[2].length) return [[w, 0]];
    return [[m[1], m[1].length / 2], [m[2], m[2].length / 2]];
  }
  function renderWord() {
    const el = $("word"); el.textContent = "";
    wordSpans = []; lastLit = -1;
    const w = cur.spell.map((s) => s[1]).join("");
    for (const [chunk, k] of groupsOf(w)) {
      const g = document.createElement("span"); g.className = "g" + (k ? " rep" : "");
      for (const ch of chunk) { const c = document.createElement("span"); c.className = "c"; c.textContent = ch; g.appendChild(c); wordSpans.push(c); }
      if (k) { const s = document.createElement("sup"); s.textContent = "×" + k; g.appendChild(s); }
      el.appendChild(g);
    }
    el.setAttribute("aria-label", `Word ${w}, ${w.length} letters`);
    const k = cur.k;
    $("formula").innerHTML = `(Ba)<sup>${k + 1}</sup> (bA)<sup>${k + 1}</sup> · <b>${4 * k + 4} letters</b>: Ba while the pair winds, bA while it unwinds`;
    $("lensCap").textContent = `winds ${k} time${k === 1 ? "" : "s"}, then unwinds`; lensTxt = "";
  }
  function spellUpdate() {
    let lit = 0;
    while (lit < cur.tauL.length && cur.tauL[lit] <= tau + 1e-9) lit++;
    if (lit === lastLit) return;
    lastLit = lit;
    wordSpans.forEach((s, i) => { s.className = "c" + (i < lit ? " on" : "") + (i === lit - 1 ? " now" : ""); });
  }
  function statusText(o) {
    if (o.status === "new") return o.extra ? "New · rung 3's word, a different orbit" : "New · not in any catalogue we could check";
    if (o.firstAtMass) return "Word known · first reported at this mass";
    return "Known · " + (o.match || "catalogued");
  }
  function select(o, keepTime) {
    cur = o;
    if (!keepTime) tau = REDUCED ? o.T * 0.62 : 0;
    $("rungTitle").innerHTML = `Rung <i>k</i> = ${o.label} <small>m₃ = ${o.m3} · T = ${fmt(o.T, 3)}</small>`;
    $("rungPill").dataset.status = o.status;
    $("rungStatus").textContent = statusText(o);
    $("legM3").textContent = String(o.m3);
    document.querySelectorAll(".rung").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.id === o.id)));
    document.querySelectorAll("#icTable tr[data-id]").forEach((r) => r.classList.toggle("sel", r.dataset.id === o.id));
    renderGhost(); renderWord(); renderHero(); renderFinger(true);
  }

  /* ladder rail */
  function thumb(cv, o) {
    fit(cv);
    const ctx = cv.getContext("2d"), d = DPR(), v = viewFor(o, cv.width, cv.height, 0.8);
    ctx.clearRect(0, 0, cv.width, cv.height);
    ghost(ctx, o, v, d, 0.75, false);
  }
  function renderRail() {
    const box = $("rungs"); box.textContent = "";
    let list = byMass(mass);
    if (!isNarrow()) list = list.slice().reverse();
    for (const o of list) {
      const b = document.createElement("button");
      b.type = "button"; b.className = "rung"; b.dataset.id = o.id; b.setAttribute("aria-pressed", String(o.id === cur.id));
      const chip = o.status === "new" ? (o.extra ? '<span class="chip extra">new · 2nd shape</span>' : '<span class="chip new">new</span>') : o.firstAtMass ? '<span class="chip">known word</span>' : '<span class="chip">catalogued</span>';
      b.innerHTML = `<canvas aria-hidden="true"></canvas><span class="meta"><span class="kk"><i>k</i>${o.label}</span><span class="tt">T = ${fmt(o.T, 2)} · ${o.spell.length} letters</span>${chip}</span>`;
      b.setAttribute("aria-label", `Rung ${o.label}, period ${fmt(o.T, 2)}, ${o.status === "new" ? "new" : "catalogued"}`);
      b.addEventListener("click", () => select(o));
      box.appendChild(b);
      requestAnimationFrame(() => thumb(b.querySelector("canvas"), o));
    }
    $("ladderNote").textContent = `m₃ = ${mass}`;
    const nNew = list.filter((o) => o.status === "new").length;
    const foot = { 0.1: "Rung 2 here came out of the m₃ = 0.1 grid search; its word is catalogued at m₃ = 0.5.", 0.15: "Rung 3 came out of the grid search at this mass; the others were followed here from m₃ = 0.25 and land on orbits the grid also found.", 0.2: "Rung 3 came out of the grid search at this mass; the others were followed here from m₃ = 0.25.", 0.3: "Rungs 3, 4 and 5 came out of the grid search at this mass; rungs 1 and 2 were followed here from m₃ = 0.25.", 0.25: "The first mass we searched. Rungs climb with an almost constant step of 9.0 in period. 3′ is a second, far less stable orbit with the rung-3 word.", 0.5: "A mass in the 2018 catalogue, which lists rungs 1 and 2 only.", 0.75: "Rung 1 here matches the catalogue to all 10 published digits. Rung 2 is not listed at this mass." }[mass] || "";
    $("ladderFoot").textContent = `${list.length} rungs, ${nNew} new. ${foot}`;
  }
  function renderMassSeg() {
    const seg = $("massSeg"); seg.textContent = "";
    for (const m of MASSES) {
      const b = document.createElement("button"); b.type = "button"; b.textContent = "m₃ " + m; b.setAttribute("aria-pressed", String(m === mass));
      b.addEventListener("click", () => {
        mass = m; seg.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
        const list = byMass(m);
        const pick = list.find((o) => o.k === cur.k && !o.extra) || list.find((o) => o.status === "new") || list[0];
        renderRail(); select(pick);
      });
      seg.appendChild(b);
    }
  }
  function goTo(o) {
    if (Math.abs(o.m3 - mass) > 1e-9) {
      mass = MASSES.find((m) => Math.abs(m - o.m3) < 1e-9) ?? mass;
      $("massSeg").querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", String(x.textContent === "m₃ " + mass)));
      renderRail();
    }
    select(o);
    document.querySelector(".plate").scrollIntoView({ behavior: REDUCED ? "auto" : "smooth", block: "start" });
  }

  /* controls */
  const playIcon = () => $("playIcon").setAttribute("d", playing ? "M4 3h3v10H4zM9 3h3v10H9z" : "M4 2.5v11l9-5.5z");
  $("play").addEventListener("click", () => { playing = !playing; playIcon(); $("play").setAttribute("aria-label", playing ? "Pause" : "Play"); });
  document.querySelectorAll("[data-speed]").forEach((b) => b.addEventListener("click", () => {
    speed = +b.dataset.speed; document.querySelectorAll("[data-speed]").forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
  }));
  playIcon();
  if (!playing) $("play").setAttribute("aria-label", "Play");

  /* ---------------- small multiples ---------------- */
  const multIds = ["L1@0.5", "L2@0.5", "L3@0.5", "L4@0.5", "L5@0.5"].filter((id) => ORB[id]);
  const mults = [];
  {
    const box = $("multiples");
    for (const id of multIds) {
      const o = ORB[id], el = document.createElement("div");
      el.className = "mult"; el.dataset.status = o.status;
      const w = o.spell.map((s) => s[1]).join("");
      el.innerHTML = `<canvas role="img" aria-label="Rung ${o.k} at m3 = 0.5"></canvas><div class="cap"><b>k = ${o.k}</b><span class="chip ${o.status === "new" ? "new" : ""}">${o.status === "new" ? "new" : "catalogued"}</span></div><div class="cap"><span class="mono">T = ${fmt(o.T, 2)}</span><span class="mono">T* = ${fmt(o.Tstar, 2)}</span></div><div class="w">(Ba)<sup>${o.k + 1}</sup>(bA)<sup>${o.k + 1}</sup> · ${o.k} turn${o.k === 1 ? "" : "s"}</div>`;
      el.querySelector("canvas").addEventListener("click", () => goTo(o));
      el.querySelector("canvas").style.cursor = "pointer";
      box.appendChild(el);
      mults.push({ o, cv: el.querySelector("canvas"), cache: null, w });
    }
    const f = $("periodFacts"), L = byMass(0.25).filter((o) => !o.extra);
    const steps = [];
    for (let i = 1; i < L.length; i++) steps.push(L[i].T - L[i - 1].T);
    const tail = steps.slice(2);
    const facts = [
      [`${fmt(Math.min(...tail), 2)}–${fmt(Math.max(...tail), 2)}`, "added to the period by each rung from k = 3 to 11 at m₃ = 0.25"],
      [`4k + 4`, "letters in the word of rung k: two more line-ups per extra loop, twice"],
      [`${L.length} rungs`, "found at m₃ = 0.25, up to a period of " + fmt(L[L.length - 1].T, 1)],
    ];
    f.innerHTML = facts.map(([n, d]) => `<div class="fact"><span class="n">${n}</span><span class="d">${d}</span></div>`).join("");
  }
  let phase = 0.0;
  function renderMults() {
    const d = DPR();
    for (const m of mults) {
      const resized = fit(m.cv);
      const w = m.cv.width, h = m.cv.height;
      if (resized || !m.cache) {
        m.cache = document.createElement("canvas"); m.cache.width = w; m.cache.height = h;
        ghost(m.cache.getContext("2d"), m.o, viewFor(m.o, w, h, 0.8), d, 0.3, false);
      }
      const ctx = m.cv.getContext("2d"), v = viewFor(m.o, w, h, 0.8);
      ctx.clearRect(0, 0, w, h); ctx.drawImage(m.cache, 0, 0);
      const t = m.o.t0 + phase * m.o.T, p = new Float64Array(6);
      syzygyFx(ctx, m.o, v, t, d * 0.7, false);
      tails(ctx, m.o, v, t, TAIL * 0.8, d, 0.6);
      posAt(m.o, t, p); bodies(ctx, v, p, m.o.m3, d, 3.6);
    }
  }

  /* ---------------- fingerprint ---------------- */
  const fOrbit = $("fOrbit"), fSphere = $("fSphere");
  let fCache = null;
  function renderFinger(force) {
    if (!fOrbit.offsetParent) return;
    const d = DPR();
    const rz = fit(fOrbit); fit(fSphere);
    const v = viewFor(cur, fOrbit.width, fOrbit.height, 0.78);
    if (rz || force || !fCache || fCache.o !== cur) {
      fCache = { o: cur, cv: document.createElement("canvas") };
      fCache.cv.width = fOrbit.width; fCache.cv.height = fOrbit.height;
      ghost(fCache.cv.getContext("2d"), cur, v, d, 0.22, true);
      $("fOrbitCap").textContent = `rung ${cur.label} · m₃ = ${cur.m3}`;
      $("fSphereCap").textContent = `${cur.syzygies.length} line-ups per period → ${cur.spell.length} letters`;
    }
    const ctx = fOrbit.getContext("2d"), t = cur.t0 + tau, p = new Float64Array(6);
    ctx.clearRect(0, 0, fOrbit.width, fOrbit.height); ctx.drawImage(fCache.cv, 0, 0);
    syzygyFx(ctx, cur, v, t, d, true);
    tails(ctx, cur, v, t, TAIL, d, 0.85);
    posAt(cur, t, p); bodies(ctx, v, p, cur.m3, d, 5);
    drawSphere(fSphere, cur, t, fYaw, fPitch, { r: 0.38, big: true });
  }
  {
    const cell = $("fSphereCell"); let px = 0, py = 0;
    cell.addEventListener("pointerdown", (e) => { dragging = true; px = e.clientX; py = e.clientY; cell.setPointerCapture(e.pointerId); cell.style.cursor = "grabbing"; });
    cell.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      fYaw += (e.clientX - px) * 0.01; fPitch = Math.max(-1.4, Math.min(1.4, fPitch + (e.clientY - py) * 0.01)); px = e.clientX; py = e.clientY;
      if (!playing) renderFinger();
    });
    const up = () => { dragging = false; cell.style.cursor = ""; };
    cell.addEventListener("pointerup", up); cell.addEventListener("pointercancel", up);
  }

  /* ---------------- family chart ---------------- */
  const SVGNS = "http://www.w3.org/2000/svg";
  function el(tag, attrs, parent, text) {
    const e = document.createElementNS(SVGNS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (text != null) e.textContent = text;
    if (parent) parent.appendChild(e);
    return e;
  }
  function familyChart() {
    const W = 960, H = 540, m = { l: 60, r: 120, t: 24, b: 52 };
    const X0 = 0.05, X1 = 1.05, Y0 = 0, Y1 = 42;
    const X = (x) => m.l + ((x - X0) / (X1 - X0)) * (W - m.l - m.r), Y = (y) => H - m.b - ((y - Y0) / (Y1 - Y0)) * (H - m.t - m.b);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Scale-free period of each rung against the mass of the third star" });
    const searched = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5].filter((q) => D.sweeps.some((s) => Math.abs(s.m3 - q) < 1e-9) || q === 0.25 || q === 0.1);
    for (const q of searched) el("rect", { x: X(q) - 6, y: m.t, width: 12, height: H - m.t - m.b, class: "band", rx: 3 }, svg);
    for (let y = 0; y <= Y1; y += 5) { el("line", { x1: m.l, x2: W - m.r, y1: Y(y), y2: Y(y), class: "grid" }, svg); el("text", { x: m.l - 10, y: Y(y) + 4, "text-anchor": "end" }, svg, y); }
    for (let x = 0.1; x <= 1.0001; x += 0.1) { el("line", { x1: X(x), x2: X(x), y1: H - m.b, y2: H - m.b + 5, class: "axis" }, svg); el("text", { x: X(x), y: H - m.b + 19, "text-anchor": "middle" }, svg, fmt(x, 1)); }
    el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: "axis" }, svg);
    el("text", { x: (m.l + W - m.r) / 2, y: H - 10, "text-anchor": "middle", class: "ttl" }, svg, "MASS OF STAR 3, m₃ (stars 1 and 2 have mass 1)");
    el("text", { x: 16, y: (m.t + H - m.b) / 2, transform: `rotate(-90 16 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle", class: "ttl" }, svg, "SCALE-FREE PERIOD T*");
    // the published catalogue: a barcode at each mass it covers
    const cat = D.catalogue.filter((r) => r[0] <= 1.05 && r[1] < Y1);
    const g0 = el("g", {}, svg);
    for (const [mm, ts] of cat) el("line", { x1: X(mm) - 7, x2: X(mm) + 7, y1: Y(ts), y2: Y(ts), stroke: C.known, "stroke-opacity": 0.28, "stroke-width": 1 }, g0);
    for (const [mm, lab] of [[0.5, "catalogued"], [0.75, "catalogued"], [1, "equal masses"]]) el("text", { x: X(mm), y: m.t - 8, "text-anchor": "middle", fill: C.muted, "font-size": 10 }, svg, lab);
    // branches
    for (const br of D.branches) {
      const neu = br.k >= 3, col = neu ? C.neu : C.known;
      const pts = br.pts.filter((p) => p[0] >= X0 && p[0] <= X1 && p[2] <= Y1);
      if (pts.length < 2) continue;
      const d = pts.map((p, i) => `${i ? "L" : "M"}${X(p[0]).toFixed(1)},${Y(p[2]).toFixed(1)}`).join("");
      if (neu) el("path", { d, fill: "none", stroke: col, "stroke-opacity": 0.18, "stroke-width": 7, "stroke-linecap": "round", "stroke-linejoin": "round" }, svg);
      el("path", { d, fill: "none", stroke: col, "stroke-width": neu ? 1.8 : 1.4, "stroke-opacity": neu ? 0.95 : 0.75, "stroke-linecap": "round", "stroke-linejoin": "round" }, svg);
      const end = pts.reduce((a, b) => (b[0] > a[0] ? b : a));
      el("text", { x: X(end[0]) + 8, y: Y(end[2]) + 4, class: neu ? "lab-new" : "", fill: neu ? C.neu : C.muted, "font-size": 12 }, svg, `k = ${br.k}`);
    }
    // rungs we computed exactly
    const tip = $("tip");
    for (const o of D.orbits) {
      const x = X(o.m3), y = Y(o.Tstar);
      if (y < m.t) continue;
      const g = el("g", { class: "pt", tabindex: 0, role: "button", "aria-label": `Rung ${o.label} at m3 = ${o.m3}, ${o.status}` }, svg);
      if (o.status === "new") {
        el("circle", { cx: x, cy: y, r: 9, fill: o.extra ? C.pink : C.neu, "fill-opacity": 0.16 }, g);
        el("circle", { cx: x, cy: y, r: 4.2, fill: o.extra ? C.pink : C.neu }, g);
      } else el("circle", { cx: x, cy: y, r: 4.2, fill: "#070b18", stroke: C.known, "stroke-width": 1.6 }, g);
      const show = (e) => { tip.hidden = false; tip.innerHTML = `rung ${o.label} · m₃ ${o.m3}<br><span class="k">T*</span> ${fmt(o.Tstar, 3)} <span class="k">T</span> ${fmt(o.T, 3)}<br><span class="k">${o.status === "new" ? "new" : o.match || "catalogued"}</span>`; tip.style.left = e.clientX + "px"; tip.style.top = e.clientY + "px"; };
      g.addEventListener("pointermove", show); g.addEventListener("pointerleave", () => (tip.hidden = true));
      g.addEventListener("click", () => { tip.hidden = true; goTo(o); });
      g.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); goTo(o); } });
    }
    const a1 = ORB["L1@1"];
    if (a1) {
      el("line", { x1: X(1) - 3, y1: Y(a1.Tstar) + 6, x2: X(1) - 14, y2: Y(a1.Tstar) + 62, stroke: C.muted, "stroke-width": 1 }, svg);
      el("text", { x: X(1) - 16, y: Y(a1.Tstar) + 76, "text-anchor": "end", fill: C.soft, "font-size": 11 }, svg, "rung 1 here is Li & Liao's I.A i.c. 3");
    }
    el("text", { x: X(0.25), y: Y(Y1) + 12, "text-anchor": "middle", fill: C.neu, "font-size": 10 }, svg, "searched");
    $("familyChart").appendChild(svg);
  }

  /* ---------------- stability chart ---------------- */
  function stabChart() {
    const W = 960, H = 380, m = { l: 60, r: 130, t: 18, b: 48 };
    const X = (k) => m.l + ((k - 0.5) / 11.5) * (W - m.l - m.r), Ly0 = 0.5, Ly1 = 3.2;
    const Y = (v) => H - m.b - ((Math.log10(v) - Ly0) / (Ly1 - Ly0)) * (H - m.t - m.b);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Largest Floquet multiplier of each rung" });
    for (const v of [3, 10, 30, 100, 300, 1000]) { el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: "grid" }, svg); el("text", { x: m.l - 10, y: Y(v) + 4, "text-anchor": "end" }, svg, v); }
    for (let k = 1; k <= 11; k++) el("text", { x: X(k), y: H - m.b + 20, "text-anchor": "middle" }, svg, k);
    el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: "axis" }, svg);
    el("text", { x: (m.l + W - m.r) / 2, y: H - 8, "text-anchor": "middle", class: "ttl" }, svg, "RUNG k");
    el("text", { x: 16, y: (m.t + H - m.b) / 2, transform: `rotate(-90 16 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle", class: "ttl" }, svg, "GROWTH PER PERIOD (LOG)");
    const series = [[0.1, "#8fd0ff"], [0.15, "#6f7da6"], [0.2, "#6f7da6"], [0.3, "#6f7da6"], [0.25, C.neu], [0.5, C.accent], [0.75, C.glow[0]]];
    for (const [ms, col] of series) {
      const pts = byMass(ms).filter((o) => !o.extra && o.maxMult).sort((a, b) => a.k - b.k);
      if (!pts.length) continue;
      if (pts.length > 1) el("path", { d: pts.map((o, i) => `${i ? "L" : "M"}${X(o.k).toFixed(1)},${Y(o.maxMult).toFixed(1)}`).join(""), fill: "none", stroke: col, "stroke-width": 1.6, "stroke-opacity": 0.85 }, svg);
      for (const o of pts) {
        const c = el("circle", { cx: X(o.k), cy: Y(o.maxMult), r: 4, fill: o.status === "new" ? col : "#070b18", stroke: col, "stroke-width": 1.6, class: "pt" }, svg);
        c.addEventListener("click", () => goTo(o));
        el("title", {}, c, `rung ${o.k}, m3 = ${o.m3}: ×${fmt(o.maxMult, 1)} per period`);
      }
      const lastP = pts[pts.length - 1];
      if (![0.15, 0.2, 0.3].includes(ms)) el("text", { x: X(lastP.k) + 10, y: Y(lastP.maxMult) + 4, fill: col, "font-size": 12 }, svg, `m₃ = ${ms}`);
    }
    const x3 = ORB["X3@0.25"];
    if (x3 && x3.maxMult) {
      const c = el("circle", { cx: X(3), cy: Y(x3.maxMult), r: 4.5, fill: C.pink, class: "pt" }, svg);
      c.addEventListener("click", () => goTo(x3));
      el("text", { x: X(3) + 10, y: Y(x3.maxMult) + 4, fill: C.pink, "font-size": 12 }, svg, `3′ (×${Math.round(x3.maxMult)})`);
    }
    $("stabChart").appendChild(svg);
  }

  /* ---------------- return-distance map ---------------- */
  const RAMP = [[0, [255, 246, 220]], [0.18, [255, 199, 106]], [0.42, [226, 80, 138]], [0.66, [90, 42, 158]], [0.88, [18, 23, 53]], [1, [7, 10, 22]]];
  function ramp(t) {
    for (let i = 1; i < RAMP.length; i++) if (t <= RAMP[i][0]) {
      const [a, ca] = RAMP[i - 1], [b, cb] = RAMP[i], f = (t - a) / (b - a);
      return [0, 1, 2].map((c) => ca[c] + f * (cb[c] - ca[c]));
    }
    return RAMP[RAMP.length - 1][1];
  }
  const MAPS = D.maps.map((M) => ({ ...M, dq: b64(M.d, Uint8Array), Tq: b64(M.T, Uint8Array) }));
  let mapIdx = 0;
  function renderMap() {
    const M = MAPS[mapIdx]; if (!M) return;
    const base = $("mapBase"), N = M.N;
    base.width = N; base.height = N;
    const ctx = base.getContext("2d"), img = ctx.createImageData(N, N);
    for (let r = 0; r < N; r++) for (let c = 0; c < N; c++) {
      const q = M.dq[r * N + c], o = ((N - 1 - r) * N + c) * 4;
      if (!q) { img.data[o] = 2; img.data[o + 1] = 4; img.data[o + 2] = 10; img.data[o + 3] = 255; continue; }
      const L = ((q - 1) / 253) * 4 - 3, t = Math.max(0, Math.min(1, (L + 2) / 2.15));
      let col = ramp(Math.pow(t, 2.1));
      if (M.Tq[r * N + c] === 255) col = col.map((x, i) => x * 0.42 + [6, 9, 22][i] * 0.5);
      img.data[o] = col[0]; img.data[o + 1] = col[1]; img.data[o + 2] = col[2]; img.data[o + 3] = 255;
    }
    ctx.putImageData(img, 0, 0);
    renderMapOver();
  }
  function mapMarkers() {
    const M = MAPS[mapIdx], out = [];
    const sw = D.sweeps.find((s) => Math.abs(s.m3 - M.m3) < 1e-9);
    if (sw) for (const o of sw.orbits) out.push({ v1: o.v1, v2: o.v2, T: o.T, kind: o.ladder ? (o.known ? "known" : "new") : o.known ? "other" : "uncat", canon: o.canon, match: o.match });
    for (const o of byMass(M.m3)) if (!out.some((q) => Math.abs(q.v1 - o.v1) < 1e-4 && Math.abs(q.v2 - o.v2) < 1e-4)) out.push({ v1: o.v1, v2: o.v2, T: o.T, kind: o.status === "new" ? "new" : "known", canon: o.canon, orbit: o });
    for (const q of out) q.orbit = q.orbit || byMass(M.m3).find((o) => Math.abs(o.v1 - q.v1) < 1e-4 && Math.abs(o.v2 - q.v2) < 1e-4);
    return out;
  }
  let hoverMarker = null;
  function renderMapOver(mx, my) {
    const over = $("mapOver"); fit(over);
    const ctx = over.getContext("2d"), w = over.width, h = over.height, d = DPR(), M = MAPS[mapIdx];
    ctx.clearRect(0, 0, w, h);
    const X = (v1) => (v1 / M.vmax) * w, Y = (v2) => h - (v2 / M.vmax) * h;
    // energy boundary
    ctx.strokeStyle = "rgba(236,235,230,0.25)"; ctx.setLineDash([3 * d, 4 * d]); ctx.lineWidth = 1 * d;
    ctx.beginPath(); ctx.arc(0, h, w, -Math.PI / 2, 0); ctx.stroke(); ctx.setLineDash([]);
    for (const q of mapMarkers()) {
      const x = X(q.v1), y = Y(q.v2), hot = hoverMarker === q;
      const col = q.kind === "new" ? [142, 240, 207] : q.kind === "uncat" ? [255, 155, 179] : [236, 235, 230];
      const r = (q.kind === "new" ? 8 : 6) * d * (hot ? 1.35 : 1);
      ctx.strokeStyle = rgba(col, q.kind === "other" ? 0.55 : 0.95); ctx.lineWidth = (q.kind === "new" ? 1.8 : 1.2) * d;
      ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.stroke();
      if (q.kind === "new") { ctx.fillStyle = rgba(col, 0.18); ctx.fill(); }
    }
    ctx.font = `500 ${Math.round(10 * d)}px "IBM Plex Mono", monospace`; ctx.fillStyle = "rgba(205,210,223,0.7)";
    ctx.fillText("v₂ ↑", 8 * d, 16 * d); ctx.textAlign = "right"; ctx.fillText("v₁ →", w - 8 * d, h - 8 * d); ctx.textAlign = "left";
    if (mx != null) { ctx.strokeStyle = "rgba(170,184,255,0.35)"; ctx.lineWidth = 1 * d; ctx.beginPath(); ctx.moveTo(mx, 0); ctx.lineTo(mx, h); ctx.moveTo(0, my); ctx.lineTo(w, my); ctx.stroke(); }
  }
  function renderMapSeg() {
    const seg = $("mapSeg"); seg.textContent = "";
    MAPS.forEach((M, i) => {
      const b = document.createElement("button"); b.type = "button"; b.textContent = `m₃ = ${M.m3}`; b.setAttribute("aria-pressed", String(i === mapIdx));
      b.addEventListener("click", () => { mapIdx = i; seg.querySelectorAll("button").forEach((x, j) => x.setAttribute("aria-pressed", String(j === i))); renderMap(); });
      seg.appendChild(b);
    });
  }
  {
    const box = $("map"), ro = $("mapReadout");
    const at = (e) => {
      const r = box.getBoundingClientRect(), M = MAPS[mapIdx], d = DPR();
      const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
      const v1 = fx * M.vmax, v2 = (1 - fy) * M.vmax;
      let best = null, bd = 14;
      for (const q of mapMarkers()) { const dd = Math.hypot((q.v1 / M.vmax - fx) * r.width, ((1 - q.v2 / M.vmax) - fy) * r.height); if (dd < bd) { bd = dd; best = q; } }
      return { v1, v2, fx, fy, best, M, d, r };
    };
    box.addEventListener("pointermove", (e) => {
      const { v1, v2, fx, fy, best, M, d, r } = at(e);
      hoverMarker = best;
      renderMapOver(fx * r.width * d, fy * r.height * d);
      const c = Math.min(M.N - 1, Math.max(0, Math.floor(fx * M.N))), row = Math.min(M.N - 1, Math.max(0, Math.floor((1 - fy) * M.N)));
      const q = M.dq[row * M.N + c], tq = M.Tq[row * M.N + c];
      const E = (v1 * v1 + v2 * v2) * (1 + 2 / M.m3) - (0.5 + 2 * M.m3);
      let html = `<span>v₁ <b>${fmt(v1, 4)}</b> · v₂ <b>${fmt(v2, 4)}</b></span><span>energy <b>${fmt(E, 4)}</b></span>`;
      if (!q) html += `<span>not bound: the stars fly apart</span>`;
      else {
        const dd = Math.pow(10, ((q - 1) / 253) * 4 - 3);
        html += `<span>closest return <b>${dd < 0.0015 ? "< 0.0015" : dd.toPrecision(2)}</b>${tq === 255 ? "" : ` at t ≈ <b>${fmt((tq / 254) * M.Tmax, 1)}</b>`}</span>`;
        if (tq === 255) html += `<span>one star was thrown out before t = ${M.Tmax}</span>`;
      }
      if (best) {
        const lab = best.orbit ? `rung ${best.orbit.label}` : "periodic orbit";
        html += `<span style="margin-top:4px;color:${best.kind === "new" ? "var(--new)" : best.kind === "uncat" ? "var(--g3)" : "var(--ink)"}">${lab} · T = ${fmt(best.T, 3)} · ${best.kind === "new" ? "new" : best.kind === "uncat" ? "word not in catalogue (unverified)" : "catalogued word"}${best.orbit ? " · click to view" : ""}</span>`;
      }
      ro.innerHTML = html;
      box.style.cursor = best && best.orbit ? "pointer" : "crosshair";
    });
    box.addEventListener("pointerleave", () => { hoverMarker = null; renderMapOver(); });
    box.addEventListener("click", (e) => { const { best } = at(e); if (best && best.orbit) goTo(best.orbit); });
  }

  /* ---------------- search numbers, verdicts, table ---------------- */
  function fillNumbers() {
    const nf = (x) => x.toLocaleString("en-US");
    const sw = D.sweeps;
    const cand = sw.reduce((s, q) => s + (q.candidates || 0), 0), jobs = sw.reduce((s, q) => s + (q.jobs || 0), 0);
    const distinct = sw.reduce((s, q) => s + q.orbits.length, 0), unc = sw.reduce((s, q) => s + q.orbits.filter((o) => !o.known).length, 0);
    const scanned = D.scan ? D.scan.total : 3518620;
    $("nScan").textContent = nf(scanned) + " orbits";
    $("stScanned").textContent = (scanned / 1e6).toFixed(1) + "M";
    const nm = sw.length, words = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"];
    $("stScanD").textContent = `starting conditions integrated across ${words[nm] || nm} masses, 800 × 800 launch velocities each`;
    $("nBranches").textContent = D.branches.length + " branches";
    $("nCand").textContent = nf(cand) + " near-returns";
    $("nRef").textContent = nf(jobs) + " refined";
    $("nDistinct").textContent = nf(distinct) + " distinct orbits";
    $("nNew").textContent = nf(unc) + " not catalogued";
    const neu = D.orbits.filter((o) => o.status === "new");
    $("stNewOrbits").textContent = neu.length;
    $("stNewWords").textContent = new Set(neu.map((o) => o.canon)).size;
    $("builtOn").textContent = "Data built " + D.built;
  }
  function fillVerdicts() {
    const V = D.verify || [
      { title: "Periodic", tag: "pending", text: "Independent re-integration with a different method." },
      { title: "Distinct", tag: "pending", text: "Independent word computation." },
      { title: "Unpublished", tag: "pending", text: "Independent catalogue comparison." },
    ];
    $("verdicts").innerHTML = V.map((v) => `<div class="verdict"><h4>${v.title}<span class="${v.warn ? "warn" : ""}">${v.tag}</span></h4><p>${v.text}</p></div>`).join("");
  }
  function compactWord(o) { return `(Ba)<sup>${o.k + 1}</sup>(bA)<sup>${o.k + 1}</sup>`; }
  function fillTable() {
    const rows = D.orbits.slice().sort((a, b) => a.m3 - b.m3 || a.k - b.k || (a.extra ? 1 : -1));
    const head = "<tr><th>Rung</th><th>m₃</th><th>Status</th><th>v₁</th><th>v₂</th><th>T</th><th>T*</th><th>Word</th><th>Growth</th><th>Closure</th><th></th></tr>";
    const body = rows.map((o) => `<tr data-id="${o.id}"><td>${o.label}</td><td>${o.m3}</td><td class="${o.status === "new" ? "ok" : ""}">${o.status === "new" ? "new" : o.firstAtMass ? "new at mass" : "catalogued"}</td><td>${o.v1.toFixed(10)}</td><td>${o.v2.toFixed(10)}</td><td>${o.T.toFixed(8)}</td><td>${o.Tstar.toFixed(4)}</td><td class="w">${compactWord(o)}</td><td>${o.maxMult ? "×" + fmt(o.maxMult, 1) : "—"}</td><td class="cl">—</td><td><button type="button" data-act="check">Check</button> <button type="button" data-act="copy">Copy</button> <button type="button" data-act="view">View</button></td></tr>`).join("");
    $("icTable").innerHTML = `<table><thead>${head}</thead><tbody>${body}</tbody></table>`;
    $("icTable").addEventListener("click", (e) => {
      const b = e.target.closest("button"); if (!b) return;
      const tr = b.closest("tr"), o = ORB[tr.dataset.id];
      if (b.dataset.act === "view") goTo(o);
      else if (b.dataset.act === "copy") {
        const s = `m3=${o.m3} v1=${o.v1} v2=${o.v2} T=${o.T}`;
        const done = () => toast("Copied the starting values for rung " + o.label);
        try { navigator.clipboard.writeText(s).then(done, () => fallbackCopy(s, done)); } catch (err) { fallbackCopy(s, done); }
      } else if (b.dataset.act === "check") {
        const cell = tr.querySelector(".cl"); cell.textContent = "integrating…"; cell.className = "cl";
        setTimeout(() => {
          const t0 = performance.now(), r = closure(o), ms = performance.now() - t0;
          cell.textContent = `${r.toExponential(1)} (${Math.round(ms)} ms)`; cell.className = "cl " + (r < 1e-6 ? "ok" : "bad");
        }, 30);
      }
    });
  }
  function fallbackCopy(s, done) {
    const ta = document.createElement("textarea"); ta.value = s; document.body.appendChild(ta); ta.select();
    try { document.execCommand("copy"); done(); } catch (e) { toast("Select the text and copy it: " + s); }
    ta.remove();
  }
  function closure(o) { // re-integrate with the Observatory solver and measure |Y(T) - Y(0)|
    const m = [1, 1, o.m3], x = [-1, 0, 0, 1, 0, 0, 0, 0, 0], v = [o.v1, o.v2, 0, o.v1, o.v2, 0, (-2 * o.v1) / o.m3, (-2 * o.v2) / o.m3, 0];
    const sim = new window.NBody(m, x, v, { tol: 1e-12 });
    while (sim.t < o.T) sim.step();
    const st = new Float64Array(18); sim.stateAt(o.T, st);
    const y0 = [...x, ...v]; let s = 0;
    for (let i = 0; i < 18; i++) s += (st[i] - y0[i]) ** 2;
    return Math.sqrt(s);
  }
  let toastTimer = 0;
  function toast(s) { const t = $("toast"); t.textContent = s; t.hidden = false; clearTimeout(toastTimer); toastTimer = setTimeout(() => (t.hidden = true), 2200); }

  /* ---------------- loop ---------------- */
  const vis = { stage: true, mult: false, finger: false };
  if ("IntersectionObserver" in window) {
    const io = new IntersectionObserver((es) => { for (const e of es) vis[e.target.dataset.vis] = e.isIntersecting; }, { rootMargin: "80px" });
    stage.dataset.vis = "stage"; $("multiples").dataset.vis = "mult"; document.querySelector(".finger").dataset.vis = "finger";
    io.observe(stage); io.observe($("multiples")); io.observe(document.querySelector(".finger"));
  } else { vis.mult = vis.finger = true; }
  let last = performance.now();
  function frame(now) {
    const dt = Math.min(0.05, (now - last) / 1000); last = now;
    if (playing) {
      tau += dt * speed * RATE;
      if (tau > cur.T) { tau -= cur.T; lastLit = -1; }
      phase = (phase + dt / 16) % 1;
      if (!dragging) fYaw += dt * 0.08;
    }
    if (vis.stage) renderHero();
    if (vis.mult && (playing || !mults[0].cache)) renderMults();
    if (vis.finger) renderFinger();
    requestAnimationFrame(frame);
  }

  function layout() { drawSky(sky); renderGhost(); renderHero(); if (MAPS.length) renderMapOver(); }
  let wasNarrow = isNarrow();
  addEventListener("resize", () => {
    clearTimeout(layout.t);
    layout.t = setTimeout(() => { layout(); if (isNarrow() !== wasNarrow) { wasNarrow = isNarrow(); renderRail(); } }, 120);
  });

  renderMassSeg(); renderRail(); fillNumbers(); fillVerdicts(); fillTable();
  familyChart(); stabChart();
  if (MAPS.length) { renderMapSeg(); renderMap(); } else $("map").hidden = true;
  drawSky(sky); select(cur);
  if (REDUCED) phase = 0.62;
  renderMults();
  requestAnimationFrame(frame);
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { renderHero(); });
})();
