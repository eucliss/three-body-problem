"""Shape-space returns: compare (r12,r13,r23, d/dt of each) at time t with the same invariants of Y0,
allowing the 1<->2 relabelling (m1=m2). A hit at t<T means the orbit is periodic modulo an isometry/relabel."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize_scalar
from nb import ic, integrate
def inv(Y):
    P = Y[:6].reshape(3, -1, order='C') if Y.ndim == 1 else None
def invariants(Yg):
    x = Yg[0:6]; v = Yg[6:12]
    out = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        dx = x[2*j] - x[2*i]; dy = x[2*j+1] - x[2*i+1]
        dvx = v[2*j] - v[2*i]; dvy = v[2*j+1] - v[2*i+1]
        r = np.hypot(dx, dy)
        out += [r, (dx * dvx + dy * dvy) / r, (dx * dvy - dy * dvx) / r]  # r, rdot, r*thetadot
    return np.array(out)
claims = json.load(open(sys.argv[1])); ids = sys.argv[2:]
for c in claims:
    if ids and c['id'] not in ids: continue
    m = (1.0, 1.0, c['m3']); Y0 = ic(c['m3'], c['v1'], c['v2']); T = c['T']
    sol = integrate(m, Y0, T)
    tg = np.linspace(0, T, 200001); I = invariants(sol.sol(tg)); I0 = invariants(Y0)
    I0s = I0[[0, 1, 2, 6, 7, 8, 3, 4, 5]]  # swap 1<->2: r12 same, r13<->r23
    I0s = I0s.copy(); I0s[2] *= -1  # angular term of pair (1,2) flips sign under relabel
    hits = []
    for name, ref in [('same', I0), ('swap12', I0s)]:
        for sgn, sr, lab in [(1, 1, ''), (-1, 1, '+reflect'), (-1, -1, '+trev'), (1, -1, '+trev+reflect')]:
            r2 = ref.copy(); r2[[2, 5, 8]] *= sgn; r2[[1, 4, 7]] *= sr
            R = np.max(np.abs(I - r2[:, None]), axis=0)
            sel = (tg > 0.01 * T) & (tg < 0.99 * T)
            j = np.argmin(np.where(sel, R, np.inf))
            hits.append(f"{name}{lab}:{R[j]:.1e}@{tg[j]/T:.4f}")
    print(f"{c['id']:<28}", ' '.join(hits), flush=True)
