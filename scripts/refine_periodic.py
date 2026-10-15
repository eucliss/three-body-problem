"""Refine the published (5-digit) Šuvakov–Dmitrašinović periodic orbits by Newton shooting.

Unknowns (p1, p2, T); residual = state(T) - state(0). Solved with Levenberg–Marquardt
on top of the regularized Bulirsch–Stoer solver. Writes data/periodic_orbits.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from threebody import Solver, presets  # noqa: E402


def residual(params):
    p1, p2, T = params
    m, x, v = presets.suvakov(p1, p2)
    tr = Solver(m, x, v, tol=1e-13).integrate(T)
    return np.concatenate([(tr.x[-1] - x)[:, :2].ravel(), (tr.v[-1] - v)[:, :2].ravel()])


def main(names):
    out_path = Path(__file__).resolve().parents[1] / "data" / "periodic_orbits.json"
    out_path.parent.mkdir(exist_ok=True)
    results = json.loads(out_path.read_text()) if out_path.exists() else {}
    for name in names:
        p1, p2, T = presets.SUVAKOV[name]
        r0 = np.abs(residual([p1, p2, T])).max()
        sol = least_squares(residual, [p1, p2, T], xtol=1e-15, ftol=1e-15, gtol=1e-15,
                            diff_step=1e-7, method="lm", max_nfev=400)
        r1 = np.abs(sol.fun).max()
        results[name] = {"p1": sol.x[0], "p2": sol.x[1], "T": sol.x[2],
                         "published": [p1, p2, T], "residual_published": r0, "residual": r1}
        print(f"{name:16s} published residual {r0:9.2e} -> refined {r1:9.2e}  "
              f"p=({sol.x[0]:.12f}, {sol.x[1]:.12f}) T={sol.x[2]:.10f}", flush=True)
        out_path.write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main(sys.argv[1:] or list(presets.SUVAKOV))
