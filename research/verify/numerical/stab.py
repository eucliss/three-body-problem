"""Monodromy matrix via variational equations (own code), Floquet multipliers,
and an independent Gauss-Newton refinement of (v1, v2, T) within the slice.
usage: python3 -I stab.py <claims.json> <out.json> [ids...]"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from nb import ic, make_rhs, energy


def make_var_rhs(m):
    f = make_rhs(m)
    pairs = [(0, 1), (0, 2), (1, 2)]

    def rhs(t, Z):
        Y = Z[:12]
        Phi = Z[12:].reshape(12, 12)
        dY = f(t, Y)
        P = Y[:6].reshape(3, 2)
        J = np.zeros((12, 12))
        J[0:6, 6:12] = np.eye(6)
        A = np.zeros((6, 6))  # d acc / d pos
        for i, j in pairs:
            d = P[j] - P[i]
            r2 = d @ d
            r = np.sqrt(r2)
            K = (np.eye(2) / r ** 3 - 3.0 * np.outer(d, d) / r ** 5)
            # a_i += m_j * d / r^3 ; d a_i / d r_j = m_j K ; d a_i/ d r_i = -m_j K
            A[2*i:2*i+2, 2*j:2*j+2] += m[j] * K
            A[2*i:2*i+2, 2*i:2*i+2] -= m[j] * K
            # a_j += m_i * (-d) / r^3 ; d a_j / d r_i = m_i K ; d a_j / d r_j = -m_i K
            A[2*j:2*j+2, 2*i:2*i+2] += m[i] * K
            A[2*j:2*j+2, 2*j:2*j+2] -= m[i] * K
        J[6:12, 0:6] = A
        return np.concatenate([dY, (J @ Phi).ravel()])
    return rhs


def monodromy(m, Y0, T, rtol=1e-13, atol=1e-13):
    Z0 = np.concatenate([Y0, np.eye(12).ravel()])
    sol = solve_ivp(make_var_rhs(m), (0, T), Z0, method='DOP853', rtol=rtol, atol=atol)
    Z = sol.y[:, -1]
    return Z[:12], Z[12:].reshape(12, 12)


def dY0(m3):
    a = np.zeros(12); a[6] = a[8] = 1.0; a[10] = -2.0 / m3
    b = np.zeros(12); b[7] = b[9] = 1.0; b[11] = -2.0 / m3
    return a, b


def run(c):
    m3 = c['m3']; m = (1.0, 1.0, m3)
    v1, v2, T = c['v1'], c['v2'], c['T']
    out = {'id': c['id']}
    t0 = time.time()
    f = make_rhs(m)
    hist = []
    for it in range(3):
        Y0 = ic(m3, v1, v2)
        YT, M = monodromy(m, Y0, T)
        F = YT - Y0
        hist.append(float(np.max(np.abs(F))))
        if it == 0:
            ev = np.linalg.eigvals(M)
            out['multipliers'] = sorted([[float(e.real), float(e.imag)] for e in ev], key=lambda z: -abs(complex(*z)))
            out['max_abs_multiplier'] = float(max(abs(ev)))
            out['det_M'] = float(np.linalg.det(M))
            # nontrivial: drop the 8 closest to 1
            evs = sorted(ev, key=lambda e: abs(e - 1))
            nt = evs[8:]
            out['nontrivial'] = [[float(e.real), float(e.imag)] for e in nt]
            out['nontrivial_abs'] = [float(abs(e)) for e in nt]
            out['linearly_stable'] = bool(all(abs(abs(e) - 1) < 1e-3 for e in nt))
            out['trivial_dist_from_1'] = [float(abs(e - 1)) for e in evs[:8]]
        a, b = dY0(m3)
        Jm = np.column_stack([M @ a - a, M @ b - b, f(0, YT)])
        dx, *_ = np.linalg.lstsq(Jm, -F, rcond=None)
        if it == 0:
            out['newton_step0'] = dx.tolist()
            out['cond_J'] = float(np.linalg.cond(Jm))
            sv = np.linalg.svd(Jm, compute_uv=False)
            out['sv_J'] = sv.tolist()
        v1 += dx[0]; v2 += dx[1]; T += dx[2]
        if hist[-1] < 2e-11:
            break
    Y0 = ic(m3, v1, v2)
    sol = solve_ivp(f, (0, T), Y0, method='DOP853', rtol=1e-13, atol=1e-13)
    hist.append(float(np.max(np.abs(sol.y[:, -1] - Y0))))
    out['newton_res_hist'] = hist
    out['refined'] = {'v1': v1, 'v2': v2, 'T': T}
    out['shift_from_claim'] = {'dv1': v1 - c['v1'], 'dv2': v2 - c['v2'], 'dT': T - c['T']}
    E = energy(Y0, m)
    out['Tstar_refined'] = T * abs(E) ** 1.5
    out['wall_s'] = time.time() - t0
    return out


if __name__ == '__main__':
    claims = json.load(open(sys.argv[1]))
    ids = sys.argv[3:]
    if ids:
        claims = [c for c in claims if c['id'] in ids]
    res = []
    for c in claims:
        o = run(c)
        res.append(o)
        s = o['shift_from_claim']
        print(f"{o['id']:<28} |lam|max={o['max_abs_multiplier']:.3e} nontriv|.|={['%.4g' % x for x in o['nontrivial_abs']]} stable={o['linearly_stable']} "
              f"trivial|e-1|max={max(o['trivial_dist_from_1']):.1e} newton res={['%.1e' % h for h in o['newton_res_hist']]} "
              f"shift=({s['dv1']:.1e},{s['dv2']:.1e},{s['dT']:.1e}) cond={o['cond_J']:.1e} t={o['wall_s']:.0f}s", flush=True)
        json.dump(res, open(sys.argv[2], 'w'), indent=1)
