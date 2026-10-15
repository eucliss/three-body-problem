/* Regularized three-body solver: Bulirsch–Stoer extrapolation of the logarithmic-Hamiltonian
 * (LogH) leapfrog (Mikkola & Tanikawa 1999; Preto & Tremaine 1999), integrated in the triangle
 * variables (the sides s_i) so close pairs keep full precision. Same algorithm as
 * threebody/integrators.py. Works in the browser (global `NBody`) and in Node (module.exports).
 * Positions/velocities are flat arrays [x0,y0,z0, x1,y1,z1, x2,y2,z2]. */
(function (root) {
  "use strict";

  const SEQ = [2, 4, 6, 8, 10, 12, 14, 16];

  /* x^(1/n) from multiplications and comparisons only. The built-in power function is not
     correctly rounded and differs between JavaScript engines; using only IEEE-754 exact
     operations (+ - * / sqrt) everywhere in the integration makes every engine take
     bit-identical steps, so a fate-map pixel computed in Node replays exactly in any browser. */
  function nthRoot(x, n) {
    if (!(x > 0)) return 0;
    x = Math.min(x, 1e12);
    let lo = 0, hi = x > 1 ? 1e4 : 1;
    for (let i = 0; i < 64; i++) {
      const mid = 0.5 * (lo + hi);
      let q = 1;
      for (let j = 0; j < n; j++) q *= mid;
      if (q > x) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
  }
  const KMAX = SEQ.length;
  const K_TARGET = 4;

  /* Three-body solver in triangle variables. The state is the three sides of the triangle,
   *   s1 = x3 - x2,  s2 = x1 - x3,  s3 = x2 - x1   (s1 + s2 + s3 = 0),
   * their velocities w_i = ds_i/dt, the barycentre R and t. Each side is integrated directly,
   * so a tight pair keeps full relative precision wherever it is: absolute coordinates would
   * only resolve a separation r near |x| ~ 0.5 to ulp(0.5)/r, which breaks the regularization
   * for r ~ 1e-9. The LogH leapfrog in these variables is the same map as in absolute
   * coordinates (drift every body by dt v, kick by (h/U) a), just evaluated on differences:
   *   d2 s_i/dt2 = -G (M - m_i) s_i/s_i^3 + G m_i sum_{j != i} s_j/s_j^3      (Lagrange 1772).
   * Absolute positions/velocities are rebuilt for output in this.y = [x(9), v(9), t]. */
  const ZL = 22; // s(9), w(9), R(3), t

  class NBody {
    constructor(m, x, v, opts) {
      opts = opts || {};
      if (m.length !== 3) throw new Error("NBody integrates exactly three bodies");
      this.N = 3;
      this.m = Float64Array.from(m);
      const M = m[0] + m[1] + m[2];
      this.M = M;
      this.G = opts.G ?? 1;
      const eps = opts.eps ?? 0;
      this.eps2 = eps * eps;
      this.tol = opts.tol ?? 1e-11;
      this.n3 = 9;
      this.len = 19; // absolute output layout: x, v, t
      // barycentre position and (constant) velocity
      this.V = new Float64Array(3);
      const z = new Float64Array(ZL);
      for (let c = 0; c < 3; c++) {
        z[18 + c] = (m[0] * x[c] + m[1] * x[3 + c] + m[2] * x[6 + c]) / M;
        this.V[c] = (m[0] * v[c] + m[1] * v[3 + c] + m[2] * v[6 + c]) / M;
        z[c] = x[6 + c] - x[3 + c]; z[3 + c] = x[c] - x[6 + c]; z[6 + c] = x[3 + c] - x[c];
        z[9 + c] = v[6 + c] - v[3 + c]; z[12 + c] = v[c] - v[6 + c]; z[15 + c] = v[3 + c] - v[c];
      }
      z[21] = opts.t ?? 0;
      this.z = z;
      this.prevZ = Float64Array.from(z);
      this.y = new Float64Array(19);
      this.toAbs(z, this.y);
      this.prev = Float64Array.from(this.y);
      this.t = z[21];
      // scratch
      this._a = new Float64Array(9);
      this._S = new Float64Array(9);
      this._W = new Float64Array(9);
      this._R = new Float64Array(3);
      this._abs = new Float64Array(19);
      this.table = [];
      for (let k = 0; k < KMAX; k++) {
        const row = [];
        for (let j = 0; j <= k; j++) row.push(new Float64Array(ZL));
        this.table.push(row);
      }
      this.steps = 0;
      this.rejects = 0;
      this.forceEvals = 0;
      this.failRun = 0;
      this.lastDt = 0;
      this.lastH = 0; this.lastMid = false;
      this.E0 = this.energy();
      this.B = -this.E0;
      // initial fictitious step ≈ 1% of the shortest pairwise dynamical time × U
      const U = this.sideAcc(z, this._a);
      const pair = [[1, 2], [2, 0], [0, 1]]; // side i joins these two bodies
      let tau = Infinity;
      for (let i = 0; i < 3; i++) {
        const mu = this.G * (m[pair[i][0]] + m[pair[i][1]]);
        if (mu <= 0) continue;
        const r = Math.sqrt(z[3 * i] * z[3 * i] + z[3 * i + 1] * z[3 * i + 1] + z[3 * i + 2] * z[3 * i + 2]);
        tau = Math.min(tau, Math.sqrt((r * r * r) / mu));
      }
      this.H = 1e-2 * U * (isFinite(tau) ? tau : 1);
    }

    /* side accelerations of the sides held in S[0..9) into a; returns U = G Σ m_i m_j / r_ij */
    sideAcc(S, a) {
      const m = this.m, G = this.G, e2 = this.eps2, M = this.M;
      const q = [0, 0, 0], u = [0, 0, 0]; // 1/r^3 and 1/r per side
      for (let i = 0; i < 3; i++) {
        const r2 = S[3 * i] * S[3 * i] + S[3 * i + 1] * S[3 * i + 1] + S[3 * i + 2] * S[3 * i + 2] + e2;
        const ir = 1 / Math.sqrt(r2);
        u[i] = ir; q[i] = ir * ir * ir;
      }
      for (let c = 0; c < 3; c++) {
        const t0 = S[c] * q[0], t1 = S[3 + c] * q[1], t2 = S[6 + c] * q[2];
        a[c] = G * (-(M - m[0]) * t0 + m[0] * (t1 + t2));
        a[3 + c] = G * (-(M - m[1]) * t1 + m[1] * (t0 + t2));
        a[6 + c] = G * (-(M - m[2]) * t2 + m[2] * (t0 + t1));
      }
      this.forceEvals++;
      return G * (m[1] * m[2] * u[0] + m[2] * m[0] * u[1] + m[0] * m[1] * u[2]);
    }

    /* kinetic energy from side velocities W[0..9): Σ_{i<j} m_i m_j |v_i - v_j|^2 / 2M + M|V|^2/2 */
    kineticW(W) {
      const m = this.m, V = this.V;
      const k1 = W[0] * W[0] + W[1] * W[1] + W[2] * W[2];
      const k2 = W[3] * W[3] + W[4] * W[4] + W[5] * W[5];
      const k3 = W[6] * W[6] + W[7] * W[7] + W[8] * W[8];
      return (m[1] * m[2] * k1 + m[2] * m[0] * k2 + m[0] * m[1] * k3) / (2 * this.M) + 0.5 * this.M * (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
    }

    potential(z) {
      return this.sideAcc(z || this.z, this._a);
    }

    energy(z) {
      z = z || this.z;
      return this.kineticW(z.subarray(9, 18)) - this.potential(z);
    }

    /* absolute [x(9), v(9), t] from triangle variables */
    toAbs(z, out) {
      const m = this.m, M = this.M, V = this.V;
      for (let c = 0; c < 3; c++) {
        const s1 = z[c], s2 = z[3 + c], s3 = z[6 + c], w1 = z[9 + c], w2 = z[12 + c], w3 = z[15 + c], R = z[18 + c];
        out[c] = R + (m[2] * s2 - m[1] * s3) / M;
        out[3 + c] = R + (m[0] * s3 - m[2] * s1) / M;
        out[6 + c] = R + (m[1] * s1 - m[0] * s2) / M;
        out[9 + c] = V[c] + (m[2] * w2 - m[1] * w3) / M;
        out[12 + c] = V[c] + (m[0] * w3 - m[2] * w1) / M;
        out[15 + c] = V[c] + (m[1] * w1 - m[0] * w2) / M;
      }
      out[18] = z[21];
      return out;
    }

    /* n substeps of the LogH leapfrog over fictitious time H: D(h/2) K(h) D(h) … K(h) D(h/2) */
    leap(z0, H, n, out) {
      const S = this._S, W = this._W, R = this._R, a = this._a, B = this.B, V = this.V;
      for (let i = 0; i < 9; i++) { S[i] = z0[i]; W[i] = z0[9 + i]; }
      for (let c = 0; c < 3; c++) R[c] = z0[18 + c];
      let t = z0[21];
      const h = H / n;
      let w = this.kineticW(W) + B;
      if (!(w > 0)) { out[0] = NaN; return; }
      let dt = (0.5 * h) / w;
      for (let i = 0; i < 9; i++) S[i] += dt * W[i];
      for (let c = 0; c < 3; c++) R[c] += dt * V[c];
      t += dt;
      for (let s = 0; s < n; s++) {
        const U = this.sideAcc(S, a);
        const f = h / U;
        for (let i = 0; i < 9; i++) W[i] += f * a[i];
        w = this.kineticW(W) + B;
        if (!(w > 0)) { out[0] = NaN; return; }
        dt = (s < n - 1 ? h : 0.5 * h) / w;
        for (let i = 0; i < 9; i++) S[i] += dt * W[i];
        for (let c = 0; c < 3; c++) R[c] += dt * V[c];
        t += dt;
      }
      for (let i = 0; i < 9; i++) { out[i] = S[i]; out[9 + i] = W[i]; }
      for (let c = 0; c < 3; c++) out[18 + c] = R[c];
      out[21] = t;
    }

    /* Gragg modified midpoint in physical time (fallback when the LogH denominator fails, and
       for landing exactly on a requested time) */
    midpoint(z0, H, n, out) {
      const a = this._a, V = this.V;
      const z1 = new Float64Array(ZL), z2 = new Float64Array(ZL), tmp = new Float64Array(ZL), d = new Float64Array(ZL);
      const f = (z, dz) => {
        this.sideAcc(z, a);
        for (let i = 0; i < 9; i++) { dz[i] = z[9 + i]; dz[9 + i] = a[i]; }
        for (let c = 0; c < 3; c++) dz[18 + c] = V[c];
        dz[21] = 1;
      };
      const h = H / n;
      z1.set(z0);
      f(z0, d);
      for (let i = 0; i < ZL; i++) z2[i] = z0[i] + h * d[i];
      for (let s = 1; s < n; s++) {
        f(z2, d);
        for (let i = 0; i < ZL; i++) { tmp[i] = z1[i] + 2 * h * d[i]; z1[i] = z2[i]; z2[i] = tmp[i]; }
      }
      f(z2, d);
      for (let i = 0; i < ZL; i++) out[i] = 0.5 * (z1[i] + z2[i] + h * d[i]);
    }

    attempt(H, scheme, z0) {
      z0 = z0 || this.z;
      const T = this.table, tol = this.tol;
      let err = Infinity;
      for (let k = 0; k < KMAX; k++) {
        const n = SEQ[k];
        scheme.call(this, z0, H, n, T[k][0]);
        if (!isFinite(T[k][0][0])) return { ok: false, k, err: Infinity };
        for (let j = 1; j <= k; j++) {
          const q = n / SEQ[k - j];
          const r = q * q - 1;
          const A = T[k][j - 1], Bp = T[k - 1][j - 1], C = T[k][j];
          for (let i = 0; i < ZL; i++) C[i] = A[i] + (A[i] - Bp[i]) / r;
        }
        if (k >= 1) {
          const Y = T[k][k], Z = T[k][k - 1];
          let s = 0;
          for (let i = 0; i < ZL; i++) {
            const sc = tol + tol * Math.max(Math.abs(z0[i]), Math.abs(Y[i]));
            const e = (Y[i] - Z[i]) / sc;
            s += e * e;
          }
          err = Math.sqrt(s / ZL);
          if (!isFinite(err)) return { ok: false, k, err: Infinity };
          if (err <= 1) return { ok: true, k, err, y: Y };
        }
      }
      return { ok: false, k: KMAX - 1, err };
    }

    static nextH(H, k, err, ok) {
      if (!ok) {
        if (!isFinite(err)) return H * 0.1;
        return H * Math.min(0.5, Math.max(0.1, 0.94 * nthRoot(0.65 / err, 2 * k + 1)));
      }
      let fac = 0.94 * nthRoot(0.65 / Math.max(err, 1e-30), 2 * k + 1);
      if (k < K_TARGET) fac *= 1.5;
      else if (k > K_TARGET) fac *= 0.7;
      return H * Math.min(4, Math.max(0.2, fac));
    }

    /* one accepted step; returns physical dt. Keeps the previous state in this.prev / prevZ.
       After 40 consecutive LogH rejections (only possible when T + B ≈ U loses all precision,
       e.g. a fully dispersed system) it takes physical-time midpoint steps until one succeeds. */
    step() {
      for (let guard = 0; guard < 200; guard++) {
        const useMid = this.failRun > 40;
        const H = useMid ? this.lastDt || 1e-6 : this.H;
        const res = this.attempt(H, useMid ? this.midpoint : this.leap);
        if (res.ok) {
          this.prevZ.set(this.z);
          this.z.set(res.y);
          this.prev.set(this.y);
          this.toAbs(this.z, this.y);
          const dt = this.z[21] - this.prevZ[21];
          this.t = this.z[21];
          this.failRun = 0;
          this.lastH = H; this.lastMid = useMid;
          if (useMid) this.lastDt = NBody.nextH(H, res.k, res.err, true);
          else { this.H = NBody.nextH(H, res.k, res.err, true); this.lastDt = dt; }
          this.steps++;
          return dt;
        }
        this.rejects++;
        if (useMid) this.lastDt = NBody.nextH(H, res.k, res.err, false);
        else { this.H = NBody.nextH(H, res.k, res.err, false); this.failRun++; }
      }
      throw new Error("solver stalled");
    }

    /* cubic Hermite interpolation of absolute positions (and velocities) at tau ∈ [prev.t, t];
       for display only (error ~1e-4 over long steps); use stateAt() for exact values */
    sample(tau, xo, vo) {
      const p = this.prev, y = this.y;
      const t0 = p[18], t1 = y[18], h = t1 - t0;
      if (!(h > 0)) {
        for (let i = 0; i < 9; i++) { xo[i] = y[i]; if (vo) vo[i] = y[9 + i]; }
        return;
      }
      const s = Math.min(1, Math.max(0, (tau - t0) / h)), s2 = s * s, s3 = s2 * s;
      const h00 = 2 * s3 - 3 * s2 + 1, h10 = s3 - 2 * s2 + s, h01 = -2 * s3 + 3 * s2, h11 = s3 - s2;
      for (let i = 0; i < 9; i++)
        xo[i] = h00 * p[i] + h10 * h * p[9 + i] + h01 * y[i] + h11 * h * y[9 + i];
      if (vo) {
        const d00 = 6 * s2 - 6 * s, d10 = 3 * s2 - 4 * s + 1, d01 = -6 * s2 + 6 * s, d11 = 3 * s2 - 2 * s;
        for (let i = 0; i < 9; i++)
          vo[i] = (d00 * p[i] + d01 * y[i]) / h + d10 * p[9 + i] + d11 * y[9 + i];
      }
    }

    /* exact absolute state [x(9), v(9)] at tau ∈ [prev.t, t], from the previous accepted state.
       The landing stays regularized: physical time after a LogH step is a smooth, increasing
       function t(H) of the fictitious step, so a bracketed secant search on H ∈ [0, lastH] hits
       tau to ~1e-14 (Illinois variant); a midpoint step of negligible length closes the gap. */
    stateAt(tau, out) {
      const z0 = this.prevZ;
      let z = Float64Array.from(z0);
      if (!this.lastMid && this.lastH > 0 && tau > z0[21]) {
        let loH = 0, loT = z0[21], hiH = this.lastH, hiT = this.z[21], side = 0;
        const best = Float64Array.from(z0);
        if (hiT >= tau) {
          for (let it = 0; it < 60; it++) {
            let H = loH + ((tau - loT) * (hiH - loH)) / (hiT - loT);
            if (!(H > loH && H < hiH)) H = 0.5 * (loH + hiH);
            const res = this.attempt(H, this.leap, z0);
            if (!res.ok) break;
            const tt = res.y[21];
            if (tt <= tau) {
              best.set(res.y); loH = H; loT = tt;
              if (side === -1) hiT = tau + 0.5 * (hiT - tau);
              side = -1;
            } else {
              hiH = H; hiT = tt;
              if (side === 1) loT = tau - 0.5 * (tau - loT);
              side = 1;
            }
            if (tt <= tau && tau - best[21] < 1e-13 * Math.max(1, Math.abs(tau))) break;
            if (hiH - loH < 1e-15 * hiH) break;
          }
          z = best;
        }
      }
      let t = z[21], h = tau - t;
      for (let guard = 0; guard < 400 && tau - t > 1e-15 * Math.max(1, Math.abs(tau)); guard++) {
        h = Math.min(h, tau - t);
        const res = this.attempt(h, this.midpoint, z);
        if (!res.ok) { h *= 0.25; continue; }
        z.set(res.y);
        t = z[21];
        h = NBody.nextH(h, res.k, res.err, true);
      }
      this.toAbs(z, this._abs);
      for (let i = 0; i < 18; i++) out[i] = this._abs[i];
      return out;
    }

    get x() { return this.y.subarray(0, 9); }
    get v() { return this.y.subarray(9, 18); }

    relEnergyError() {
      const E = this.energy();
      return this.E0 !== 0 ? Math.abs((E - this.E0) / this.E0) : Math.abs(E);
    }

    angularMomentum(y) {
      y = y || this.y;
      const L = [0, 0, 0];
      for (let i = 0; i < 3; i++) {
        const k = 3 * i, m = this.m[i];
        const x = y[k], yy = y[k + 1], z = y[k + 2], vx = y[9 + k], vy = y[9 + k + 1], vz = y[9 + k + 2];
        L[0] += m * (yy * vz - z * vy); L[1] += m * (z * vx - x * vz); L[2] += m * (x * vy - yy * vx);
      }
      return L;
    }
  }

  /* ---- analysis helpers (3 bodies) ---- */

  function pairElements(x, v, m, i, j, G) {
    G = G ?? 1;
    const rx = x[3 * j] - x[3 * i], ry = x[3 * j + 1] - x[3 * i + 1], rz = x[3 * j + 2] - x[3 * i + 2];
    const wx = v[3 * j] - v[3 * i], wy = v[3 * j + 1] - v[3 * i + 1], wz = v[3 * j + 2] - v[3 * i + 2];
    const mu = G * (m[i] + m[j]);
    const r = Math.sqrt(rx * rx + ry * ry + rz * rz);
    const eps = 0.5 * (wx * wx + wy * wy + wz * wz) - mu / r;
    const hx = ry * wz - rz * wy, hy = rz * wx - rx * wz, hz = rx * wy - ry * wx;
    const h2 = hx * hx + hy * hy + hz * hz;
    const e = Math.sqrt(Math.max(0, 1 + (2 * eps * h2) / (mu * mu)));
    const a = eps < 0 ? -mu / (2 * eps) : Infinity;
    const hn = Math.sqrt(h2);
    const inc = hn > 0 ? (Math.acos(Math.max(-1, Math.min(1, hz / hn))) * 180) / Math.PI : 0;
    return { a, e, eps, r, inc };
  }

  /* Dynamical state of a triple: 'escape' | 'hierarchical' | 'interplay' | 'unbound' */
  function classify(x, v, m, G, rMin) {
    G = G ?? 1;
    let best = null;
    const pairs = [[0, 1, 2], [1, 2, 0], [0, 2, 1]];
    for (const [i, j, k] of pairs) {
      const el = pairElements(x, v, m, i, j, G);
      if (el.eps < 0 && (!best || el.a < best.el.a)) best = { i, j, k, el };
    }
    if (!best) return { state: "unbound" };
    const { i, j, k, el } = best;
    const mb = m[i] + m[j];
    const R = [0, 0, 0], V = [0, 0, 0];
    for (let c = 0; c < 3; c++) {
      R[c] = x[3 * k + c] - (m[i] * x[3 * i + c] + m[j] * x[3 * j + c]) / mb;
      V[c] = v[3 * k + c] - (m[i] * v[3 * i + c] + m[j] * v[3 * j + c]) / mb;
    }
    const Rn = Math.sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    const Ek = 0.5 * (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]) - (G * (mb + m[k])) / Rn;
    const rv = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
    const out = { pair: [i, j], third: k, a: el.a, e: el.e, inc: el.inc, R: Rn, Ek };
    if (Ek > 0 && rv > 0 && Rn > Math.max(5 * el.a, rMin || 0)) {
      out.state = "escape";
      out.vInf = Math.sqrt(2 * Ek);
    } else if (Ek < 0 && Rn > 3 * el.a) out.state = "hierarchical";
    else out.state = "interplay";
    return out;
  }

  /* Montgomery's shape sphere: unit vector; poles = Lagrange triangles, equator = collinear */
  function shapeSphere(x, m) {
    const [m1, m2, m3] = m;
    const r1 = [0, 0, 0], r2 = [0, 0, 0];
    for (let c = 0; c < 3; c++) {
      r1[c] = x[3 + c] - x[c];
      r2[c] = x[6 + c] - (m1 * x[c] + m2 * x[3 + c]) / (m1 + m2);
    }
    const mu1 = Math.sqrt((m1 * m2) / (m1 + m2)), mu2 = Math.sqrt(((m1 + m2) * m3) / (m1 + m2 + m3));
    const z1 = r1.map((q) => q * mu1), z2 = r2.map((q) => q * mu2);
    const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const cr = [z1[1] * z2[2] - z1[2] * z2[1], z1[2] * z2[0] - z1[0] * z2[2], z1[0] * z2[1] - z1[1] * z2[0]];
    const w1 = dot(z1, z1) - dot(z2, z2), w2 = 2 * dot(z1, z2);
    // signed area: orient by +z when the triangle is (nearly) planar in xy, else by its own normal
    const crn = Math.hypot(cr[0], cr[1], cr[2]);
    const w3 = 2 * (Math.abs(cr[2]) > 1e-9 * crn ? Math.sign(cr[2]) : 1) * crn;
    const n = Math.hypot(w1, w2, w3);
    return n > 0 ? [w1 / n, w2 / n, w3 / n] : [0, 0, 0];
  }

  const api = { NBody, classify, pairElements, shapeSphere };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else Object.assign(root, api);
})(typeof window !== "undefined" ? window : globalThis);
