"""Look for times t in (0,T) where Y(t) equals a symmetry image of Y(0) (relabel 1<->2, reflections, time reversal)."""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize_scalar
from nb import ic, integrate
def images(Y):
    out = {}
    P = Y[:6].reshape(3, 2); V = Y[6:].reshape(3, 2)
    for sw, rx, ry, tr in itertools.product([0, 1], [1, -1], [1, -1], [1, -1]):
        if (sw, rx, ry, tr) == (0, 1, 1, 1): continue
        Q = P.copy(); W = V.copy()
        if sw: Q = Q[[1, 0, 2]]; W = W[[1, 0, 2]]
        Q = Q * np.array([rx, ry]); W = W * np.array([rx, ry]) * tr
        out[f"sw{sw}rx{rx}ry{ry}tr{tr}"] = np.concatenate([Q.ravel(), W.ravel()])
    return out
claims = json.load(open(sys.argv[1])); ids = sys.argv[2:]
for c in claims:
    if ids and c['id'] not in ids: continue
    m = (1.0, 1.0, c['m3']); Y0 = ic(c['m3'], c['v1'], c['v2']); T = c['T']
    sol = integrate(m, Y0, T)
    tg = np.linspace(0, T, 200001); Yg = sol.sol(tg)
    hits = []
    for k, S in images(Y0).items():
        R = np.max(np.abs(Yg - S[:, None]), axis=0)
        sel = (tg > 1e-3 * T) & (tg < (1 - 1e-3) * T)
        j = np.argmin(np.where(sel, R, np.inf))
        if R[j] < 1e-2:
            g = lambda t: np.max(np.abs(sol.sol(t) - S))
            r = minimize_scalar(g, bounds=(tg[j - 2], tg[j + 2]), method='bounded', options={'xatol': 1e-13})
            hits.append(f"{k}@t/T={r.x / T:.6f}:{r.fun:.1e}")
    print(f"{c['id']:<28}", hits if hits else 'no symmetry-image return', flush=True)
