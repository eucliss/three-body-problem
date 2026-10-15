"""Integrators for the gravitational N-body equations.

Main solver: ``Solver`` — Bulirsch–Stoer extrapolation of the time-symmetric
logarithmic-Hamiltonian (LogH) leapfrog of Mikkola & Tanikawa (1999) and
Preto & Tremaine (1999). Integrating in the fictitious time s with

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

from .physics import accel_potential, energy, kinetic

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
    x: np.ndarray  # (steps, N, 3)
    v: np.ndarray
    stats: dict = field(default_factory=dict)


class Solver:
    """Regularized high-precision N-body solver (LogH leapfrog + Bulirsch–Stoer)."""

    def __init__(self, m, x, v, G=1.0, eps=0.0, tol=1e-12, t=0.0):
        self.m = np.asarray(m, float)
        self.N = len(self.m)
        self.G = float(G)
        self.eps2 = float(eps) ** 2
        self.rtol = self.atol = float(tol)
        self.x = np.array(x, float).reshape(self.N, 3)
        self.v = np.array(v, float).reshape(self.N, 3)
        self.t = float(t)
        self.E0 = energy(self.x, self.v, self.m, self.G, self.eps2)
        self.B = -self.E0
        self.steps = 0
        self.rejects = 0
        self.force_evals = 0
        # initial fictitious step ≈ 1% of the shortest pairwise dynamical time, times U
        _, U = accel_potential(self.x, self.m, self.G, self.eps2)
        tau = np.inf
        for i in range(self.N):
            for j in range(i + 1, self.N):
                r = np.linalg.norm(self.x[i] - self.x[j])
                mu = self.G * (self.m[i] + self.m[j])
                if mu > 0:
                    tau = min(tau, np.sqrt(r**3 / mu))
        self.H = 1e-2 * U * (tau if np.isfinite(tau) else 1.0)

    # --- packing -----------------------------------------------------------
    def _pack(self, x, v, t):
        return np.concatenate([x.ravel(), v.ravel(), [t]])

    def _unpack(self, y):
        n3 = 3 * self.N
        return y[:n3].reshape(self.N, 3), y[n3 : 2 * n3].reshape(self.N, 3), y[-1]

    @property
    def y(self):
        return self._pack(self.x, self.v, self.t)

    # --- schemes -----------------------------------------------------------
    def _logh_leapfrog(self, y0, H, n):
        """n substeps of the LogH leapfrog D(h/2) K(h) D(h) ... K(h) D(h/2)."""
        m, G, e2, B = self.m, self.G, self.eps2, self.B
        x, v, t = self._unpack(y0)
        x = x.copy()
        v = v.copy()
        h = H / n
        dt = 0.5 * h / (kinetic(v, m) + B)
        x += dt * v
        t += dt
        for i in range(n):
            a, U = accel_potential(x, m, G, e2)
            v += (h / U) * a
            dt = (h if i < n - 1 else 0.5 * h) / (kinetic(v, m) + B)
            x += dt * v
            t += dt
        self.force_evals += n
        return self._pack(x, v, t)

    def _midpoint(self, y0, H, n):
        """Gragg's modified midpoint rule in physical time (unregularized)."""
        m, G, e2 = self.m, self.G, self.eps2
        n3 = 3 * self.N

        def f(y):
            x = y[:n3].reshape(self.N, 3)
            a, _ = accel_potential(x, m, G, e2)
            return np.concatenate([y[n3 : 2 * n3], a.ravel(), [1.0]])

        h = H / n
        z0 = y0
        z1 = y0 + h * f(y0)
        for _ in range(n - 1):
            z0, z1 = z1, z0 + 2 * h * f(z1)
        self.force_evals += n + 1
        return 0.5 * (z0 + z1 + h * f(z1))

    # --- stepping ----------------------------------------------------------
    def step(self, H_max=None):
        """Take one accepted regularized step. Returns the physical time advanced."""
        while True:
            H = self.H if H_max is None else min(self.H, H_max)
            y_new, k, err = _bs_attempt(self.y, H, self._logh_leapfrog, self.rtol, self.atol)
            if y_new is not None:
                x, v, t = self._unpack(y_new)
                dt = t - self.t
                self.x, self.v, self.t = x.copy(), v.copy(), float(t)
                self.H = _next_step(H, k, err, True)
                self.steps += 1
                return dt
            self.rejects += 1
            self.H = _next_step(H, k, err, False)
            if self.H < 1e-300:
                raise FloatingPointError("step size underflow")

    def _land(self, t_out):
        """State at exactly t_out (≥ current t) without disturbing the main integration:
        adaptive Bulirsch–Stoer on the modified midpoint rule, from the current state."""
        y = self.y.copy()
        t = self.t
        h = t_out - t
        while t_out - t > 0:
            h = min(h, t_out - t)
            y_new, k, err = _bs_attempt(y, h, self._midpoint, self.rtol, self.atol)
            if y_new is None:
                h *= 0.25
                continue
            y = y_new
            t += h
            h = _next_step(h, k, err, True)
        x, v, _ = self._unpack(y)
        return x, v

    def integrate(self, t_end, dt_out=None, stop=None, max_steps=10_000_000):
        """Integrate to t_end.

        dt_out : if given, record states at exact multiples of dt_out (uniform sampling);
                 otherwise record every accepted step.
        stop   : optional callback stop(t, x, v) -> bool, checked after each step.
        """
        ts, xs, vs = [self.t], [self.x.copy()], [self.v.copy()]
        next_out = self.t + dt_out if dt_out else None
        stopped = False
        while self.t < t_end and self.steps < max_steps:
            t_prev, x_prev, v_prev = self.t, self.x.copy(), self.v.copy()
            H_save = self.H
            self.step()
            if self.t > t_end:  # overshoot: land exactly on t_end from the previous state
                self.t, self.x, self.v, self.H = t_prev, x_prev, v_prev, H_save
                if dt_out:
                    while next_out < t_end - 1e-15:
                        xo, vo = self._land(next_out)
                        ts.append(next_out); xs.append(xo); vs.append(vo)
                        next_out += dt_out
                x_end, v_end = self._land(t_end)
                self.x, self.v, self.t = x_end, v_end, float(t_end)
                ts.append(self.t); xs.append(self.x.copy()); vs.append(self.v.copy())
                break
            if dt_out:
                if next_out <= self.t:
                    # outputs inside this step: land on each one from the step start
                    t_now, x_now, v_now = self.t, self.x, self.v
                    self.t, self.x, self.v = t_prev, x_prev, v_prev
                    while next_out <= t_now:
                        xo, vo = self._land(next_out)
                        ts.append(next_out); xs.append(xo); vs.append(vo)
                        next_out += dt_out
                    self.t, self.x, self.v = t_now, x_now, v_now
            else:
                ts.append(self.t); xs.append(self.x.copy()); vs.append(self.v.copy())
            if stop is not None and stop(self.t, self.x, self.v):
                stopped = True
                break
        E = energy(self.x, self.v, self.m, self.G, self.eps2)
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
