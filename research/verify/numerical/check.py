"""Per-claim numerical check. Usage: python3 -I check.py <claims.json> <out.json> [ids...]"""
import sys, json, time
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize_scalar
from nb import ic, integrate, energy, angmom, momentum, pair_dists, make_rhs


def analyze(c):
    m3 = c['m3']; m = (1.0, 1.0, m3)
    Y0 = ic(m3, c['v1'], c['v2'])
    T = c['T']
    out = {'id': c['id']}
    E0 = energy(Y0, m)
    out['E_ic'] = E0
    out['E_claim_diff'] = E0 - c['E']
    out['Tstar_recomputed'] = T * abs(E0) ** 1.5
    out['Tstar_claim_diff'] = out['Tstar_recomputed'] - c['Tstar']
    out['L0'] = angmom(Y0, m)
    out['P0'] = momentum(Y0, m).tolist()
    t0 = time.time()
    sol = integrate(m, Y0, T)
    out['wall_s'] = time.time() - t0
    out['nsteps'] = len(sol.t)
    out['status_ivp'] = sol.status
    YT = sol.y[:, -1]
    d = YT - Y0
    out['res_max'] = float(np.max(np.abs(d)))
    out['res_pos'] = float(np.max(np.abs(d[:6])))
    out['res_vel'] = float(np.max(np.abs(d[6:])))
    out['res_vec'] = d.tolist()
    # conserved quantities along the step points
    Es = np.array([energy(sol.y[:, k], m) for k in range(0, sol.y.shape[1], max(1, sol.y.shape[1] // 2000))] + [energy(YT, m)])
    out['E_drift_max'] = float(np.max(np.abs(Es - E0)))
    out['E_drift_end'] = float(energy(YT, m) - E0)
    Ls = np.array([angmom(sol.y[:, k], m) for k in range(0, sol.y.shape[1], max(1, sol.y.shape[1] // 2000))])
    out['L_max'] = float(np.max(np.abs(Ls)))
    # min pairwise distance: step points + dense grid, then local refine
    tg = np.unique(np.concatenate([sol.t, np.linspace(0, T, 400001)]))
    Yg = sol.sol(tg)
    D = pair_dists(Yg[:6].T)
    k = np.argmin(D.min(axis=1)); pair = np.argmin(D[k])
    lo = tg[max(k - 2, 0)]; hi = tg[min(k + 2, len(tg) - 1)]
    f = lambda t: pair_dists(sol.sol(t)[:6])[pair]
    r = minimize_scalar(f, bounds=(lo, hi), method='bounded', options={'xatol': 1e-12})
    out['rmin'] = float(min(r.fun, D[k, pair]))
    out['rmin_pair'] = ['12', '13', '23'][pair]
    out['rmin_t'] = float(r.x)
    # per-pair minima
    out['rmin_by_pair'] = D.min(axis=0).tolist()
    # Distance from initial state along orbit
    R = np.max(np.abs(Yg - Y0[:, None]), axis=0)
    # earlier returns at T/n
    subs = {}
    for n in range(2, 7):
        tc = T / n
        w = 0.05 * T / n
        sel = (tg > tc - w) & (tg < tc + w)
        j = np.argmin(np.where(sel, R, np.inf))
        g = lambda t: np.max(np.abs(sol.sol(t) - Y0))
        rr = minimize_scalar(g, bounds=(tg[max(j - 3, 0)], tg[min(j + 3, len(tg) - 1)]), method='bounded', options={'xatol': 1e-12})
        subs[str(n)] = {'at_T_over_n': float(g(tc)), 'min_near': float(min(rr.fun, R[j])), 't_min': float(rr.x)}
    out['sub_returns'] = subs
    # global earliest close approach to initial state in (0.02T, 0.98T)
    sel = (tg > 0.1 * T) & (tg < 0.9 * T)
    j = np.argmin(np.where(sel, R, np.inf))
    out['min_return_interior'] = float(R[j]); out['min_return_interior_t_over_T'] = float(tg[j] / T)
    # residual minimum in time near T (is T slightly off?)
    g = lambda t: np.linalg.norm(sol.sol(t) - Y0)
    out['res_norm2_at_T'] = float(np.linalg.norm(d))
    # tighter / looser tolerance comparison
    t0 = time.time()
    sol12 = integrate(m, Y0, T, rtol=1e-12, atol=1e-12, dense=False)
    out['res_max_tol12'] = float(np.max(np.abs(sol12.y[:, -1] - Y0)))
    out['YT_diff_13_vs_12'] = float(np.max(np.abs(sol12.y[:, -1] - YT)))
    sol14 = integrate(m, Y0, T, rtol=2.3e-14, atol=1e-15, dense=False)
    out['res_max_tol14'] = float(np.max(np.abs(sol14.y[:, -1] - Y0)))
    out['YT_diff_13_vs_14'] = float(np.max(np.abs(sol14.y[:, -1] - YT)))
    out['wall_total_s'] = time.time() - t0 + out['wall_s']
    return out


if __name__ == '__main__':
    claims = json.load(open(sys.argv[1]))
    ids = sys.argv[3:]
    if ids:
        claims = [c for c in claims if c['id'] in ids]
    res = []
    for c in claims:
        o = analyze(c)
        res.append(o)
        print(f"{o['id']:<30} res={o['res_max']:.2e} (tol12 {o['res_max_tol12']:.2e}, tol14 {o['res_max_tol14']:.2e}) "
              f"rmin={o['rmin']:.3e}[{o['rmin_pair']}] dE={o['E_drift_max']:.1e} L={o['L_max']:.1e} "
              f"sub=" + ",".join(f"{n}:{v['min_near']:.2e}" for n, v in o['sub_returns'].items()) +
              f" intmin={o['min_return_interior']:.2e}@{o['min_return_interior_t_over_T']:.3f} steps={o['nsteps']} t={o['wall_total_s']:.0f}s",
              flush=True)
        json.dump(res, open(sys.argv[2], 'w'), indent=1)
