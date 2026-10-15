/* Regularized N-body solver: Bulirsch–Stoer extrapolation of the logarithmic-Hamiltonian
 * (LogH) leapfrog (Mikkola & Tanikawa 1999; Preto & Tremaine 1999).
 * Same algorithm as threebody/integrators.py. Works in the browser (global `NBody`) and in
 * Node (module.exports). Positions/velocities are flat arrays [x0,y0,z0, x1,y1,z1, ...]. */
(function (root) {
  "use strict";

  const SEQ = [2, 4, 6, 8, 10, 12, 14, 16];
  const KMAX = SEQ.length;
  const K_TARGET = 4;

  class NBody {
    constructor(m, x, v, opts) {
      opts = opts || {};
      this.N = m.length;
      this.m = Float64Array.from(m);
      this.G = opts.G ?? 1;
      this.eps2 = (opts.eps ?? 0) ** 2;
      this.tol = opts.tol ?? 1e-11;
      this.t = opts.t ?? 0;
      const n3 = 3 * this.N;
      this.n3 = n3;
      this.len = 2 * n3 + 1; // packed state: x, v, t
      this.y = new Float64Array(this.len);
      this.y.set(x, 0);
      this.y.set(v, n3);
      this.y[2 * n3] = this.t;
      // scratch
      this._a = new Float64Array(n3);
      this._lx = new Float64Array(n3);
      this._lv = new Float64Array(n3);
      this.table = [];
      for (let k = 0; k < KMAX; k++) {
        const row = [];
        for (let j = 0; j <= k; j++) row.push(new Float64Array(this.len));
        this.table.push(row);
      }
      // previous accepted state, for Hermite interpolation between steps
      this.prev = new Float64Array(this.len);
      this.prev.set(this.y);
      this.E0 = this.energy();
      this.B = -this.E0;
      this.steps = 0;
      this.rejects = 0;
      this.forceEvals = 0;
      this.failRun = 0;
      this.lastDt = 0;
      // initial fictitious step ≈ 1% of the shortest pairwise dynamical time × U
      const U = this.accel(this.y, this._a);
      let tau = Infinity;
      for (let i = 0; i < this.N; i++)
        for (let j = i + 1; j < this.N; j++) {
          const mu = this.G * (this.m[i] + this.m[j]);
          if (mu <= 0) continue;
          const dx = x[3 * j] - x[3 * i], dy = x[3 * j + 1] - x[3 * i + 1], dz = x[3 * j + 2] - x[3 * i + 2];
          const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
          tau = Math.min(tau, Math.sqrt((r * r * r) / mu));
        }
      this.H = 1e-2 * U * (isFinite(tau) ? tau : 1);
    }

    /* accelerations of the positions held in p[0..3N) into a; returns U = G Σ m_i m_j / r_ij */
    accel(p, a) {
      const N = this.N, m = this.m, G = this.G, e2 = this.eps2;
      a.fill(0);
      let U = 0;
      for (let i = 0; i < N; i++) {
        const ix = 3 * i;
        for (let j = i + 1; j < N; j++) {
          const jx = 3 * j;
          const dx = p[jx] - p[ix], dy = p[jx + 1] - p[ix + 1], dz = p[jx + 2] - p[ix + 2];
          const r2 = dx * dx + dy * dy + dz * dz + e2;
          const ir = 1 / Math.sqrt(r2);
          const ir3 = ir * ir * ir;
          const si = G * m[j] * ir3, sj = G * m[i] * ir3;
          a[ix] += si * dx; a[ix + 1] += si * dy; a[ix + 2] += si * dz;
          a[jx] -= sj * dx; a[jx + 1] -= sj * dy; a[jx + 2] -= sj * dz;
          U += m[i] * m[j] * ir;
        }
      }
      this.forceEvals++;
      return G * U;
    }

    kinetic(vel, off) {
      off = off || 0;
      let T = 0;
      for (let i = 0; i < this.N; i++) {
        const k = off + 3 * i;
        T += this.m[i] * (vel[k] * vel[k] + vel[k + 1] * vel[k + 1] + vel[k + 2] * vel[k + 2]);
      }
      return 0.5 * T;
    }

    potential(p) {
      return this.accel(p || this.y, this._a);
    }

    energy(y) {
      y = y || this.y;
      return this.kinetic(y, this.n3) - this.potential(y);
    }

    /* n substeps of the LogH leapfrog over fictitious time H: D(h/2) K(h) D(h) … K(h) D(h/2) */
    leap(y0, H, n, out) {
      const n3 = this.n3, x = this._lx, v = this._lv, a = this._a, B = this.B;
      for (let i = 0; i < n3; i++) { x[i] = y0[i]; v[i] = y0[n3 + i]; }
      let t = y0[2 * n3];
      const h = H / n;
      let w = this.kinetic(v) + B;
      if (!(w > 0)) { out[0] = NaN; return; }
      let dt = (0.5 * h) / w;
      for (let i = 0; i < n3; i++) x[i] += dt * v[i];
      t += dt;
      for (let s = 0; s < n; s++) {
        const U = this.accel(x, a);
        const f = h / U;
        for (let i = 0; i < n3; i++) v[i] += f * a[i];
        w = this.kinetic(v) + B;
        if (!(w > 0)) { out[0] = NaN; return; }
        dt = (s < n - 1 ? h : 0.5 * h) / w;
        for (let i = 0; i < n3; i++) x[i] += dt * v[i];
        t += dt;
      }
      for (let i = 0; i < n3; i++) { out[i] = x[i]; out[n3 + i] = v[i]; }
      out[2 * n3] = t;
    }

    /* Gragg modified midpoint in physical time (fallback when the LogH denominator fails) */
    midpoint(y0, H, n, out) {
      const n3 = this.n3, len = this.len, a = this._a;
      const z0 = new Float64Array(len), z1 = new Float64Array(len), tmp = new Float64Array(len);
      const f = (y, d) => {
        this.accel(y, a);
        for (let i = 0; i < n3; i++) { d[i] = y[n3 + i]; d[n3 + i] = a[i]; }
        d[2 * n3] = 1;
      };
      const h = H / n, d = new Float64Array(len);
      z0.set(y0);
      f(y0, d);
      for (let i = 0; i < len; i++) z1[i] = y0[i] + h * d[i];
      for (let s = 1; s < n; s++) {
        f(z1, d);
        for (let i = 0; i < len; i++) { tmp[i] = z0[i] + 2 * h * d[i]; z0[i] = z1[i]; z1[i] = tmp[i]; }
      }
      f(z1, d);
      for (let i = 0; i < len; i++) out[i] = 0.5 * (z0[i] + z1[i] + h * d[i]);
    }

    attempt(H, scheme, y0) {
      y0 = y0 || this.y;
      const T = this.table, len = this.len, tol = this.tol;
      let err = Infinity;
      for (let k = 0; k < KMAX; k++) {
        const n = SEQ[k];
        scheme.call(this, y0, H, n, T[k][0]);
        if (!isFinite(T[k][0][0])) return { ok: false, k, err: Infinity };
        for (let j = 1; j <= k; j++) {
          const q = n / SEQ[k - j];
          const r = q * q - 1;
          const A = T[k][j - 1], Bp = T[k - 1][j - 1], C = T[k][j];
          for (let i = 0; i < len; i++) C[i] = A[i] + (A[i] - Bp[i]) / r;
        }
        if (k >= 1) {
          const Y = T[k][k], Z = T[k][k - 1];
          let s = 0;
          for (let i = 0; i < len; i++) {
            const sc = tol + tol * Math.max(Math.abs(y0[i]), Math.abs(Y[i]));
            const e = (Y[i] - Z[i]) / sc;
            s += e * e;
          }
          err = Math.sqrt(s / len);
          if (!isFinite(err)) return { ok: false, k, err: Infinity };
          if (err <= 1) return { ok: true, k, err, y: Y };
        }
      }
      return { ok: false, k: KMAX - 1, err };
    }

    static nextH(H, k, err, ok) {
      if (!ok) {
        if (!isFinite(err)) return H * 0.1;
        return H * Math.min(0.5, Math.max(0.1, 0.94 * Math.pow(0.65 / err, 1 / (2 * k + 1))));
      }
      let fac = 0.94 * Math.pow(0.65 / Math.max(err, 1e-30), 1 / (2 * k + 1));
      if (k < K_TARGET) fac *= 1.5;
      else if (k > K_TARGET) fac *= 0.7;
      return H * Math.min(4, Math.max(0.2, fac));
    }

    /* one accepted step; returns physical dt. Keeps the previous state in this.prev. */
    step() {
      for (let guard = 0; guard < 200; guard++) {
        const useMid = this.failRun > 40;
        const H = useMid ? this.lastDt || 1e-6 : this.H;
        const res = this.attempt(H, useMid ? this.midpoint : this.leap);
        if (res.ok) {
          this.prev.set(this.y);
          this.y.set(res.y);
          const dt = this.y[2 * this.n3] - this.prev[2 * this.n3];
          this.t = this.y[2 * this.n3];
          if (useMid) { this.failRun = 0; this.lastDt = NBody.nextH(H, res.k, res.err, true); }
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

    /* cubic Hermite interpolation of positions (and velocities) at time tau ∈ [prev.t, t] */
    sample(tau, xo, vo) {
      const n3 = this.n3, p = this.prev, y = this.y;
      const t0 = p[2 * n3], t1 = y[2 * n3], h = t1 - t0;
      if (!(h > 0)) {
        for (let i = 0; i < n3; i++) { xo[i] = y[i]; if (vo) vo[i] = y[n3 + i]; }
        return;
      }
      const s = Math.min(1, Math.max(0, (tau - t0) / h)), s2 = s * s, s3 = s2 * s;
      const h00 = 2 * s3 - 3 * s2 + 1, h10 = s3 - 2 * s2 + s, h01 = -2 * s3 + 3 * s2, h11 = s3 - s2;
      for (let i = 0; i < n3; i++)
        xo[i] = h00 * p[i] + h10 * h * p[n3 + i] + h01 * y[i] + h11 * h * y[n3 + i];
      if (vo) {
        const d00 = 6 * s2 - 6 * s, d10 = 3 * s2 - 4 * s + 1, d01 = -6 * s2 + 6 * s, d11 = 3 * s2 - 2 * s;
        for (let i = 0; i < n3; i++)
          vo[i] = (d00 * p[i] + d01 * y[i]) / h + d10 * p[n3 + i] + d11 * y[n3 + i];
      }
    }

    /* exact state at time tau ∈ [prev.t, t]: adaptive Bulirsch–Stoer (midpoint, physical time)
       from the previous accepted state; accurate to the solver tolerance. */
    stateAt(tau, out) {
      const n3 = this.n3;
      const y = Float64Array.from(this.prev);
      let t = y[2 * n3], h = tau - t;
      for (let guard = 0; guard < 400 && tau - t > 1e-15 * Math.max(1, Math.abs(tau)); guard++) {
        h = Math.min(h, tau - t);
        const res = this.attempt(h, this.midpoint, y);
        if (!res.ok) { h *= 0.25; continue; }
        y.set(res.y);
        t += h;
        h = NBody.nextH(h, res.k, res.err, true);
      }
      out.set(y.subarray(0, 2 * n3));
      return out;
    }

    get x() { return this.y.subarray(0, this.n3); }
    get v() { return this.y.subarray(this.n3, 2 * this.n3); }

    relEnergyError() {
      const E = this.energy();
      return this.E0 !== 0 ? Math.abs((E - this.E0) / this.E0) : Math.abs(E);
    }

    angularMomentum(y) {
      y = y || this.y;
      const n3 = this.n3, L = [0, 0, 0];
      for (let i = 0; i < this.N; i++) {
        const k = 3 * i, m = this.m[i];
        const x = y[k], yy = y[k + 1], z = y[k + 2], vx = y[n3 + k], vy = y[n3 + k + 1], vz = y[n3 + k + 2];
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
    const Rn = Math.hypot(R[0], R[1], R[2]);
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
