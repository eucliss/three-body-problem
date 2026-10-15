"""Integrators for the gravitational N-body equations.

Main solver: ``Solver`` — Bulirsch–Stoer extrapolation of the time-symmetric
logarithmic-Hamiltonian (LogH) leapfrog of Mikkola & Tanikawa (1999) and
Preto & Tremaine (1999), run in the triangle variables (the sides of the triangle the three
bodies form) so that close pairs keep full precision. Integrating in the fictitious time s with

    dt/ds = 1 / (T + B),        B = -E0  (binding energy, constant)
    dv/ds = a(x) / U(x)

makes the leapfrog *algorithmically regularized*: a two-body Kepler orbit is
reproduced exactly (only its phase drifts), even through a collision. Richardson
extrapolation of the leapfrog in (h/n)² then gives arbitrarily high order, with
adaptive step and order control. Tolerance 1e-12 typically conserves energy to
~1e-13 relative.

Also provided for comparison: classic Bulirsch–Stoer on the modified midpoint
rule (in physical time), fixed-step leapfrog (kick-drift-kick) and RK4.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .physics import accel_potential

SEQ = (2, 4, 6, 8, 10, 12, 14, 16)
K_TARGET = 4  # aim to converge on column 4 (n = 10), i.e. order ~10


def _bs_attempt(y0, H, scheme, rtol, atol):
    """One Bulirsch–Stoer macro-step: run `scheme` with n = 2, 4, 6, ... substeps and
    extrapolate to h → 0 (Aitken–Neville on h²). Returns (y, k, err) or (None, k, err)."""
    table = []
    err = np.inf
    for k, n in enumerate(SEQ):
        row = [scheme(y0, H, n)]
        for j in range(1, k + 1):
            r = (n / SEQ[k - j]) ** 2 - 1.0
            row.append(row[j - 1] + (row[j - 1] - table[k - 1][j - 1]) / r)
        table.append(row)
        if k >= 1:
            sc = atol + rtol * np.maximum(np.abs(y0), np.abs(row[k]))
            err = float(np.sqrt(np.mean(((row[k] - row[k - 1]) / sc) ** 2)))
            if not np.isfinite(err):
                return None, k, np.inf
            if err <= 1.0:
                return row[k], k, err
    return None, len(SEQ) - 1, err


def _next_step(H, k, err, accepted):
    if not accepted:
        if not np.isfinite(err):
            return H * 0.1
        return H * min(0.5, max(0.1, 0.94 * (0.65 / err) ** (1.0 / (2 * k + 1))))
    fac = 0.94 * (0.65 / max(err, 1e-30)) ** (1.0 / (2 * k + 1))
    if k < K_TARGET:
        fac *= 1.5
    elif k > K_TARGET:
        fac *= 0.7
    return H * min(4.0, max(0.2, fac))


@dataclass
class Trajectory:
    t: np.ndarray
    x: np.ndarray  # (steps, 3, 3)
    v: np.ndarray
    stats: dict = field(default_factory=dict)


class Solver:
    """Regularized high-precision three-body solver (LogH leapfrog + Bulirsch–Stoer).

    Integrates in triangle variables: the sides s1 = x3 - x2, s2 = x1 - x3, s3 = x2 - x1 and
    their velocities, plus the barycentre. Each side is a state variable of its own, so a tight
    pair keeps full relative precision wherever it is (absolute coordinates resolve a separation
    r near |x| ~ 0.5 only to ulp(0.5)/r). Side accelerations follow Lagrange's symmetric form

        d2 s_i/dt2 = -G (M - m_i) s_i/s_i^3 + G m_i sum_{j != i} s_j/s_j^3.

    The LogH leapfrog is the same map as in absolute coordinates, evaluated on differences.
    State vector z = [s (9), w (9), R (3), t].
    """

    PAIRS = ((1, 2), (2, 0), (0, 1))  # side i joins these two bodies

    def __init__(self, m, x, v, G=1.0, eps=0.0, tol=1e-12, t=0.0):
        self.m = np.asarray(m, float)
        if len(self.m) != 3:
            raise ValueError("Solver integrates exactly three bodies")
        self.N = 3
        self.M = float(self.m.sum())
        self.G = float(G)
        self.eps2 = float(eps) ** 2
        self.rtol = self.atol = float(tol)
        x = np.array(x, float).reshape(3, 3)
        v = np.array(v, float).reshape(3, 3)
        self.V = (self.m[:, None] * v).sum(0) / self.M
        R = (self.m[:, None] * x).sum(0) / self.M
        S = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
        W = np.array([v[2] - v[1], v[0] - v[2], v[1] - v[0]])
        self.z = np.concatenate([S.ravel(), W.ravel(), R, [float(t)]])
        self.steps = 0
        self.rejects = 0
        self.force_evals = 0
        self.fail_run = 0
        self.last_dt = 0.0
        self.last_H, self.last_mid = 0.0, False
        self.E0 = self.energy()
        self.B = -self.E0
        # initial fictitious step ≈ 1% of the shortest pairwise dynamical time, times U
        _, U = self._side_acc(S)
        tau = np.inf
        for i, (a, b) in enumerate(self.PAIRS):
            mu = self.G * (self.m[a] + self.m[b])
            if mu > 0:
                tau = min(tau, np.sqrt(np.linalg.norm(S[i]) ** 3 / mu))
        self.H = 1e-2 * U * (tau if np.isfinite(tau) else 1.0)

    # --- state ---------------------------------------------------------------
    @property
    def t(self):
        return float(self.z[-1])

    def _abs(self, z):
        S = z[:9].reshape(3, 3)
        W = z[9:18].reshape(3, 3)
        R = z[18:21]
        m, M = self.m, self.M
        x = np.array([R + (m[2] * S[1] - m[1] * S[2]) / M,
                      R + (m[0] * S[2] - m[2] * S[0]) / M,
                      R + (m[1] * S[0] - m[0] * S[1]) / M])
        v = np.array([self.V + (m[2] * W[1] - m[1] * W[2]) / M,
                      self.V + (m[0] * W[2] - m[2] * W[0]) / M,
                      self.V + (m[1] * W[0] - m[0] * W[1]) / M])
        return x, v

    @property
    def x(self):
        return self._abs(self.z)[0]

    @property
    def v(self):
        return self._abs(self.z)[1]

    # --- physics in triangle variables ----------------------------------------
    def _side_acc(self, S):
        """Side accelerations (3, 3) and U for sides S (3, 3)."""
        m, M, G = self.m, self.M, self.G
        r2 = np.einsum("ij,ij->i", S, S) + self.eps2
        ir = 1.0 / np.sqrt(r2)
        T = S * (ir**3)[:, None]  # s_i / s_i^3
        # sum over the OTHER two sides, added directly: total - own would cancel catastrophically
        # for a tight side, whose s/s^3 term can be ~1e18 while the tidal sum is O(1)
        others = np.array([T[1] + T[2], T[0] + T[2], T[0] + T[1]])
        a = G * (-(M - m)[:, None] * T + m[:, None] * others)
        U = G * (m[1] * m[2] * ir[0] + m[2] * m[0] * ir[1] + m[0] * m[1] * ir[2])
        self.force_evals += 1
        return a, U

    def _kinetic(self, W):
        m = self.m
        k = np.einsum("ij,ij->i", W, W)
        return (m[1] * m[2] * k[0] + m[2] * m[0] * k[1] + m[0] * m[1] * k[2]) / (2 * self.M) + 0.5 * self.M * (self.V @ self.V)

    def energy(self, z=None):
        z = self.z if z is None else z
        _, U = self._side_acc(z[:9].reshape(3, 3))
        return self._kinetic(z[9:18].reshape(3, 3)) - U

    # --- schemes -------------------------------------------------------------
    def _logh_leapfrog(self, z0, H, n):
        """n substeps of the LogH leapfrog D(h/2) K(h) D(h) ... K(h) D(h/2)."""
        S = z0[:9].reshape(3, 3).copy()
        W = z0[9:18].reshape(3, 3).copy()
        R = z0[18:21].copy()
        t = z0[21]
        B, V = self.B, self.V
        h = H / n
        w = self._kinetic(W) + B
        if not w > 0:
            return np.full_like(z0, np.nan)
        dt = 0.5 * h / w
        S += dt * W
        R += dt * V
        t += dt
        for i in range(n):
            a, U = self._side_acc(S)
            W += (h / U) * a
            w = self._kinetic(W) + B
            if not w > 0:
                return np.full_like(z0, np.nan)
            dt = (h if i < n - 1 else 0.5 * h) / w
            S += dt * W
            R += dt * V
            t += dt
        return np.concatenate([S.ravel(), W.ravel(), R, [t]])

    def _midpoint(self, z0, H, n):
        """Gragg's modified midpoint rule in physical time (unregularized)."""
        V = self.V

        def f(z):
            a, _ = self._side_acc(z[:9].reshape(3, 3))
            return np.concatenate([z[9:18], a.ravel(), V, [1.0]])

        h = H / n
        z1 = z0
        z2 = z0 + h * f(z0)
        for _ in range(n - 1):
            z1, z2 = z2, z1 + 2 * h * f(z2)
        return 0.5 * (z1 + z2 + h * f(z2))

    # --- stepping ------------------------------------------------------------
    def step(self):
        """Take one accepted regularized step. Returns the physical time advanced.
        After 40 consecutive LogH rejections (only when T + B ~ U has lost all precision, e.g. a
        fully dispersed system) it takes physical-time midpoint steps until one succeeds."""
        for _ in range(200):
            use_mid = self.fail_run > 40
            H = (self.last_dt or 1e-6) if use_mid else self.H
            z_new, k, err = _bs_attempt(self.z, H, self._midpoint if use_mid else self._logh_leapfrog, self.rtol, self.atol)
            if z_new is not None:
                dt = z_new[21] - self.z[21]
                self.z = z_new.copy()
                self.fail_run = 0
                self.last_H, self.last_mid = H, use_mid
                if use_mid:
                    self.last_dt = _next_step(H, k, err, True)
                else:
                    self.H = _next_step(H, k, err, True)
                    self.last_dt = dt
                self.steps += 1
                return dt
            self.rejects += 1
            if use_mid:
                self.last_dt = _next_step(H, k, err, False)
            else:
                self.H = _next_step(H, k, err, False)
                self.fail_run += 1
        raise FloatingPointError("solver stalled")

    def _land(self, z_start, t_out, H_hi=None, mid=False):
        """Triangle state at exactly t_out, starting from z_start, without disturbing the main
        integration. The landing stays regularized: physical time after a LogH step is a smooth,
        increasing function t(H) of the fictitious step H, so a bracketed secant search on H in
        [0, H_hi] (H_hi = the accepted step that passed t_out) hits t_out to ~1e-14; a final
        midpoint step of negligible length closes the gap exactly."""
        z = z_start
        if H_hi is not None and not mid and t_out > z_start[21]:
            lo_H, lo_t = 0.0, z_start[21]
            hi_H, hi_t = H_hi, None
            z_hi, _, _ = _bs_attempt(z_start, hi_H, self._logh_leapfrog, self.rtol, self.atol)
            if z_hi is not None and z_hi[21] >= t_out:
                hi_t = z_hi[21]
                best = z_start
                side = 0
                for _ in range(60):
                    H = lo_H + (t_out - lo_t) * (hi_H - lo_H) / (hi_t - lo_t)
                    if not lo_H < H < hi_H:
                        H = 0.5 * (lo_H + hi_H)
                    z_try, _, _ = _bs_attempt(z_start, H, self._logh_leapfrog, self.rtol, self.atol)
                    if z_try is None:
                        break
                    tt = z_try[21]
                    if tt <= t_out:
                        best, lo_H, lo_t = z_try, H, tt
                        if side == -1:  # Illinois modification against one-sided stalls
                            hi_t = t_out + 0.5 * (hi_t - t_out)
                        side = -1
                    else:
                        hi_H, hi_t = H, tt
                        if side == 1:
                            lo_t = t_out - 0.5 * (t_out - lo_t)
                        side = 1
                    if t_out - best[21] < 1e-13 * max(1.0, abs(t_out)) and tt <= t_out:
                        break
                    if hi_H - lo_H < 1e-15 * hi_H:
                        break
                z = best
        # close any remaining gap with physical-time midpoint steps (tiny when the search worked)
        z = z.copy()
        h = t_out - z[21]
        while t_out - z[21] > 1e-15 * max(1.0, abs(t_out)):
            h = min(h, t_out - z[21])
            z_new, k, err = _bs_attempt(z, h, self._midpoint, self.rtol, self.atol)
            if z_new is None:
                h *= 0.25
                continue
            z = z_new
            h = _next_step(h, k, err, True)
        return z

    def integrate(self, t_end, dt_out=None, stop=None, max_steps=10_000_000):
        """Integrate to t_end.

        dt_out : if given, record states at t0 + k*dt_out (uniform sampling, computed as
                 multiples so the grid does not drift); otherwise record every accepted step.
        stop   : optional callback stop(t, x, v) -> bool, checked after each step.
        """
        t0 = self.t
        x0, v0 = self._abs(self.z)
        ts, xs, vs = [t0], [x0], [v0]
        k_out = 1
        guard = 1e-9 * dt_out if dt_out else 0.0
        stopped = False

        def emit(z):
            xa, va = self._abs(z)
            ts.append(float(z[21])); xs.append(xa); vs.append(va)

        while self.t < t_end and self.steps < max_steps:
            z_prev, H_prev = self.z.copy(), self.H
            self.step()
            land = lambda tt: self._land(z_prev, tt, self.last_H, self.last_mid)
            if self.t > t_end:  # overshoot: land exactly on t_end from the previous state
                if dt_out:
                    while t0 + k_out * dt_out < t_end - guard:
                        emit(land(t0 + k_out * dt_out))
                        k_out += 1
                z_end = land(t_end)
                self.z, self.H = z_end, H_prev
                emit(self.z)
                break
            if dt_out:
                while t0 + k_out * dt_out <= self.t + guard:
                    emit(land(t0 + k_out * dt_out))
                    k_out += 1
            else:
                emit(self.z)
            if stop is not None and stop(self.t, *self._abs(self.z)):
                stopped = True
                break
        E = self.energy()
        stats = {
            "steps": self.steps,
            "rejects": self.rejects,
            "force_evals": self.force_evals,
            "dE_rel": abs((E - self.E0) / self.E0) if self.E0 != 0 else abs(E - self.E0),
            "stopped": stopped,
        }
        return Trajectory(np.array(ts), np.array(xs), np.array(vs), stats)


# --- simple fixed-step methods, for comparison only ------------------------------

def leapfrog(m, x, v, dt, n_steps, G=1.0, eps=0.0):
    """Kick-drift-kick leapfrog (2nd order, symplectic)."""
    m = np.asarray(m, float)
    x = np.array(x, float)
    v = np.array(v, float)
    e2 = eps * eps
    a, _ = accel_potential(x, m, G, e2)
    out = [x.copy()]
    for _ in range(n_steps):
        v += 0.5 * dt * a
        x += dt * v
        a, _ = accel_potential(x, m, G, e2)
        v += 0.5 * dt * a
        out.append(x.copy())
    return np.array(out), x, v


def rk4(m, x, v, dt, n_steps, G=1.0, eps=0.0):
    """Classic 4th-order Runge–Kutta (not symplectic: energy drifts)."""
    m = np.asarray(m, float)
    x = np.array(x, float)
    v = np.array(v, float)
    e2 = eps * eps

    def acc(p):
        return accel_potential(p, m, G, e2)[0]

    out = [x.copy()]
    for _ in range(n_steps):
        k1x, k1v = v, acc(x)
        k2x, k2v = v + 0.5 * dt * k1v, acc(x + 0.5 * dt * k1x)
        k3x, k3v = v + 0.5 * dt * k2v, acc(x + 0.5 * dt * k2x)
        k4x, k4v = v + dt * k3v, acc(x + dt * k3x)
        x = x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        v = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        out.append(x.copy())
    return np.array(out), x, v
