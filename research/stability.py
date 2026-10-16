"""Linear stability of a planar periodic three-body orbit from its 12x12 monodromy matrix.

The full monodromy has eight eigenvalues equal to 1 (translations, momenta, rotation, angular
momentum, time shift, energy) that sit in Jordan blocks, so finite-difference errors of size e
split them by ~sqrt(e) and fake an instability. We remove them exactly: restrict to perturbations
that conserve centre of mass, momentum, energy and angular momentum (a 6-dimensional subspace C),
then quotient out the two neutral directions inside C (time shift Ydot and rotation). The
remaining 4x4 block carries the nontrivial Floquet multipliers. The orbit is linearly stable when
all four lie on the unit circle.

Input JSON (stdin or file): {m: [3], y0: [12] (x1,y1,x2,y2,x3,y3,vx1,...), T, M: 12x12, ...}
"""
import json
import sys

import numpy as np


def reduced_multipliers(m, y0, M, G=1.0):
    m = np.asarray(m, float)
    y0 = np.asarray(y0, float)
    M = np.asarray(M, float)
    x = y0[:6].reshape(3, 2)
    v = y0[6:].reshape(3, 2)
    # accelerations, energy and angular-momentum gradients at y0
    a = np.zeros((3, 2))
    dU = np.zeros((3, 2))  # gradient of the potential energy -U (U>0) wrt positions
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            d = x[j] - x[i]
            r3 = np.linalg.norm(d) ** 3
            a[i] += G * m[j] * d / r3
            dU[i] -= G * m[i] * m[j] * d / r3  # d(-U)/dx_i
    ydot = np.concatenate([v.ravel(), a.ravel()])
    rot = np.concatenate([np.column_stack([-x[:, 1], x[:, 0]]).ravel(), np.column_stack([-v[:, 1], v[:, 0]]).ravel()])
    gE = np.concatenate([dU.ravel(), (m[:, None] * v).ravel()])
    gL = np.concatenate([(m[:, None] * np.column_stack([v[:, 1], -v[:, 0]])).ravel(),
                         (m[:, None] * np.column_stack([-x[:, 1], x[:, 0]])).ravel()])
    gCx = np.concatenate([np.repeat([[1, 0]], 3, 0) * m[:, None], np.zeros((3, 2))]).ravel()
    gCy = np.concatenate([np.repeat([[0, 1]], 3, 0) * m[:, None], np.zeros((3, 2))]).ravel()
    gPx = np.concatenate([np.zeros((3, 2)), np.repeat([[1, 0]], 3, 0) * m[:, None]]).ravel()
    gPy = np.concatenate([np.zeros((3, 2)), np.repeat([[0, 1]], 3, 0) * m[:, None]]).ravel()
    Gm = np.vstack([gCx, gCy, gPx, gPy, gE, gL])
    # basis of C = null space of the six conserved-quantity gradients
    _, s, Vt = np.linalg.svd(Gm)
    C = Vt[6:].T  # 12 x 6
    # neutral directions inside C, then an orthonormal complement of them within C
    N = np.column_stack([ydot, rot])
    Nc = C @ np.linalg.lstsq(C, N, rcond=None)[0]  # project into C (they already lie there)
    Q, _ = np.linalg.qr(np.column_stack([Nc, C]))
    Bn, Bu = Q[:, :2], Q[:, 2:6]
    B = np.column_stack([Bn, Bu])
    coeff = np.linalg.lstsq(B, M @ Bu, rcond=None)[0]  # 6 x 4
    Mred = coeff[2:, :]  # action on the quotient C / span(neutral)
    leak = float(np.linalg.norm(Gm @ (M @ Bu)) / max(1.0, np.linalg.norm(M @ Bu)))
    return np.linalg.eigvals(Mred), leak


def analyse(it):
    ev, leak = reduced_multipliers(it["m"], it["y0"], it["M"])
    mods = np.abs(ev)
    lam = float(mods.max())
    T = it.get("T")
    return {**{k: v for k, v in it.items() if k not in ("M", "y0")},
            "multipliers": [[float(z.real), float(z.imag)] for z in ev],
            "max_multiplier": lam,
            "stable": bool(lam < 1 + 1e-5),
            "lyapunov_per_time": float(np.log(lam)) / T if T else None,
            "conservation_leak": leak}


if __name__ == "__main__":
    src = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
    data = json.load(src)
    print(json.dumps([analyse(d) for d in data] if isinstance(data, list) else analyse(data)))
