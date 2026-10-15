"""Outcome classification and chaos diagnostics."""
from __future__ import annotations

import numpy as np

from .physics import two_body_elements


def classify(x, v, m, G=1.0, r_escape=None):
    """Decide the current dynamical state of a triple.

    Returns a dict with 'state' in {'escape', 'hierarchical', 'interplay', 'unbound'} and,
    when a binary exists, its members, semi-major axis and eccentricity. A body counts as
    escaping when it is unbound from the barycentre of the tightest bound pair, moving away
    from it, and farther than r_escape (default: 5 × the binary's semi-major axis)."""
    best = None
    for i, j, k in ((0, 1, 2), (1, 2, 0), (0, 2, 1)):
        el = two_body_elements(x[i], v[i], m[i], x[j], v[j], m[j], G)
        if el["eps"] < 0 and (best is None or el["a"] < best[3]["a"]):
            best = (i, j, k, el)
    if best is None:
        return {"state": "unbound"}
    i, j, k, el = best
    mb = m[i] + m[j]
    xb = (m[i] * x[i] + m[j] * x[j]) / mb
    vb = (m[i] * v[i] + m[j] * v[j]) / mb
    R = x[k] - xb
    V = v[k] - vb
    Rn = np.linalg.norm(R)
    Ek = 0.5 * V @ V - G * (mb + m[k]) / Rn  # specific energy of k about the binary
    limit = r_escape if r_escape is not None else 5 * el["a"]
    out = {"pair": (i, j), "third": k, "a": el["a"], "e": el["e"], "R": Rn, "Ek": Ek}
    if Ek > 0 and R @ V > 0 and Rn > limit:
        out["state"] = "escape"
        out["v_inf"] = np.sqrt(2 * Ek)
    elif Ek < 0 and Rn > 3 * el["a"]:
        out["state"] = "hierarchical"
    else:
        out["state"] = "interplay"
    return out


def phase_distance(xa, va, xb, vb):
    return float(np.sqrt(np.sum((xa - xb) ** 2) + np.sum((va - vb) ** 2)))


def lyapunov_fit(t, d, d_sat=1e-2):
    """Slope of ln(separation) vs time before the separation saturates: an estimate of the
    largest Lyapunov exponent."""
    t = np.asarray(t)
    d = np.asarray(d)
    mask = (d > 0) & (d < d_sat)
    if mask.sum() < 5:
        return np.nan
    p = np.polyfit(t[mask], np.log(d[mask]), 1)
    return p[0]
