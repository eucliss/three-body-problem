"""Independent shape-sphere / free-group-word reader for planar 3-body orbits.

Written from scratch for verification; does not use the project's code.
State y = [x1,y1,x2,y2,x3,y3, vx1,vy1,vx2,vy2,vx3,vy3], G = 1.
"""
import numpy as np
from scipy.integrate import solve_ivp

SQ2, SQ6 = np.sqrt(2.0), np.sqrt(6.0)


def rhs_factory(m):
    m1, m2, m3 = m

    def rhs(t, y):
        r1, r2, r3 = y[0:2], y[2:4], y[4:6]
        d12 = r2 - r1
        d13 = r3 - r1
        d23 = r3 - r2
        i12 = 1.0 / np.dot(d12, d12) ** 1.5
        i13 = 1.0 / np.dot(d13, d13) ** 1.5
        i23 = 1.0 / np.dot(d23, d23) ** 1.5
        a1 = m2 * d12 * i12 + m3 * d13 * i13
        a2 = -m1 * d12 * i12 + m3 * d23 * i23
        a3 = -m1 * d13 * i13 - m2 * d23 * i23
        return np.concatenate([y[6:12], a1, a2, a3])

    return rhs


def area(t, y):
    """twice the signed area of triangle 1-2-3 (proportional to +n_z)."""
    r1, r2, r3 = y[0:2], y[2:4], y[4:6]
    a = r2 - r1
    b = r3 - r1
    return a[0] * b[1] - a[1] * b[0]


def darea(y):
    r1, r2, r3 = y[0:2], y[2:4], y[4:6]
    v1, v2, v3 = y[6:8], y[8:10], y[10:12]
    a, b = r2 - r1, r3 - r1
    da, db = v2 - v1, v3 - v1
    return (da[0] * b[1] - da[1] * b[0]) + (a[0] * db[1] - a[1] * db[0])


def initial_state(m3, v1, v2):
    return np.array([-1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                     v1, v2, v1, v2, -2 * v1 / m3, -2 * v2 / m3])


def energy(y, m):
    m1, m2, m3 = m
    v = y[6:12].reshape(3, 2)
    r = y[0:6].reshape(3, 2)
    K = 0.5 * (m1 * v[0] @ v[0] + m2 * v[1] @ v[1] + m3 * v[2] @ v[2])
    U = -(m1 * m2 / np.linalg.norm(r[0] - r[1]) + m1 * m3 / np.linalg.norm(r[0] - r[2])
          + m2 * m3 / np.linalg.norm(r[1] - r[2]))
    return K + U


def shape_n(y):
    """Li & Liao 2017 shape-sphere vector n (works on (...,12) arrays)."""
    y = np.asarray(y)
    r1, r2, r3 = y[..., 0:2], y[..., 2:4], y[..., 4:6]
    rho = (r1 - r2) / SQ2
    lam = (r1 + r2 - 2 * r3) / SQ6
    R2 = (rho ** 2).sum(-1) + (lam ** 2).sum(-1)
    nx = 2 * (rho * lam).sum(-1) / R2
    ny = ((lam ** 2).sum(-1) - (rho ** 2).sum(-1)) / R2
    nz = 2 * (rho[..., 0] * lam[..., 1] - rho[..., 1] * lam[..., 0]) / R2
    return np.stack([nx, ny, nz], -1)


# punctures on the sphere for this (equal-mass Jacobi) parametrisation
PUNCT = {
    (1, 2): np.array([0.0, 1.0, 0.0]),
    (1, 3): np.array([-np.sqrt(3) / 2, -0.5, 0.0]),
    (2, 3): np.array([np.sqrt(3) / 2, -0.5, 0.0]),
}


def middle_body(y):
    r = y[0:6].reshape(3, 2)
    # direction of the (nearly) collinear line: longest side
    pairs = [(0, 1), (0, 2), (1, 2)]
    i, j = max(pairs, key=lambda p: np.linalg.norm(r[p[0]] - r[p[1]]))
    d = r[j] - r[i]
    d /= np.linalg.norm(d)
    s = r @ d
    order = np.argsort(s)
    return int(order[1]) + 1  # 1-based


def integrate(m3, v1, v2, T, rtol=1e-13, atol=1e-13, tmax=None, method="DOP853"):
    m = (1.0, 1.0, m3)
    y0 = initial_state(m3, v1, v2)
    ev = area
    ev.direction = 0
    if tmax is None:
        tmax = 1.3 * T
    sol = solve_ivp(rhs_factory(m), (0.0, tmax), y0, method=method, rtol=rtol, atol=atol,
                    events=[ev], dense_output=True)
    if sol.status != 0:
        raise RuntimeError(sol.message)
    return sol


def crossings(sol, T, tskip=1e-9):
    """Return t_s and list of (t, middle_body, dir) in window (t_s, t_s+T]."""
    te = sol.t_events[0]
    ye = sol.y_events[0]
    te_pos = te[te > tskip]
    if len(te_pos) == 0:
        raise RuntimeError("no crossing")
    t1 = te_pos[0]
    ts = 0.5 * t1
    out = []
    for t, y in zip(te, ye):
        if ts < t <= ts + T:
            d = darea(y)
            out.append((t, middle_body(y), 1 if d > 0 else -1, abs(d)))
    return ts, out


def free_reduce_cyclic(w):
    st = []
    for c in w:
        if st and st[-1] == c.swapcase():
            st.pop()
        else:
            st.append(c)
    # cyclic reduction
    while len(st) >= 2 and st[0] == st[-1].swapcase():
        st = st[1:-1]
    return "".join(st)


def inverse(w):
    return "".join(c.swapcase() for c in reversed(w))


def canon(w):
    """canonical rep. up to cyclic rotation and inversion."""
    w = free_reduce_cyclic(w)
    if not w:
        return ""
    cands = []
    for x in (w, inverse(w)):
        for k in range(len(x)):
            cands.append(x[k:] + x[:k])
    return min(cands)


def expand(word):
    """expand notation like (ab)2ABAba(BAB) , b2(ABab)2A2 ... into plain letters."""
    out = []
    i = 0
    s = word.replace(" ", "")

    def parse(i):
        res = []
        while i < len(s):
            c = s[i]
            if c == "(":
                sub, i = parse(i + 1)
                # now s[i] == ')'
                i += 1
                n, i = num(i)
                res.extend(sub * n)
            elif c == ")":
                return res, i
            elif c.isalpha():
                i += 1
                n, i = num(i)
                res.extend([c] * n)
            else:
                raise ValueError(word)
        return res, i

    def num(i):
        j = i
        while j < len(s) and s[j].isdigit():
            j += 1
        return (int(s[i:j]) if j > i else 1), j

    res, _ = parse(0)
    return "".join(res)


# ---- convention enumeration (method A: combinatorial, body-in-middle) ----
CONVENTIONS = []
for pole in [(1, 2), (1, 3), (2, 3)]:
    for abody in pole:
        for low in (+1, -1):
            CONVENTIONS.append((pole, abody, low))


def word_A(cross, conv):
    pole, abody, low = conv
    w = []
    for (t, mid, d, _) in cross:
        if mid not in pole:
            continue
        c = "a" if mid == abody else "b"
        w.append(c if d == low else c.upper())
    return "".join(w)


def word_B(sol, ts, T, conv, nsamp=400000):
    """Geometric method: real stereographic projection from the pole puncture,
    classify crossings of the image of the equator by X position."""
    pole, abody, low = conv
    P = PUNCT[pole]
    others = [k for k in PUNCT if k != pole]
    e1 = np.cross([0, 0, 1.0], P)
    tt = np.linspace(ts, ts + T, nsamp)
    Y = sol.sol(tt).T
    n = shape_n(Y)
    den = 1 - n @ P
    X = (n @ e1) / den
    Z = n[:, 2] / den
    Xp = {k: (PUNCT[k] @ e1) / (1 - PUNCT[k] @ P) for k in others}
    # which puncture is 'a': the puncture adjacent (along equator arc) to the
    # arc with abody in the middle. arc E_i joins C_ij(pole) and C_ik.
    k3 = ({1, 2, 3} - set(pole)).pop()
    apunct = tuple(sorted((abody, k3)))
    s = np.sign(Z)
    idx = np.nonzero(s[1:] != s[:-1])[0]
    w = []
    xa = Xp[apunct]
    bp = [k for k in others if k != apunct][0]
    xb = Xp[bp]
    for i in idx:
        f = Z[i] / (Z[i] - Z[i + 1])
        x = X[i] + f * (X[i + 1] - X[i])
        d = 1 if Z[i + 1] > Z[i] else -1
        # ray from a-puncture outward (away from the other puncture)
        if (x - xa) * (xa - xb) > 0:
            c = "a"
        elif (x - xb) * (xb - xa) > 0:
            c = "b"
        else:
            continue
        w.append(c if d == low else c.upper())
    return "".join(w)


def min_puncture_dist(sol, ts, T, nsamp=200000):
    tt = np.linspace(ts, ts + T, nsamp)
    n = shape_n(sol.sol(tt).T)
    return {k: float(np.min(np.linalg.norm(n - P, axis=1))) for k, P in PUNCT.items()}


def word_B_events(cross_states, conv):
    """Geometric classification at root-found crossing states (stereographic X)."""
    pole, abody, low = conv
    P = PUNCT[pole]
    others = [k for k in PUNCT if k != pole]
    e1 = np.cross([0, 0, 1.0], P)
    k3 = ({1, 2, 3} - set(pole)).pop()
    apunct = tuple(sorted((abody, k3)))
    bp = [k for k in others if k != apunct][0]
    xa = (PUNCT[apunct] @ e1) / (1 - PUNCT[apunct] @ P)
    xb = (PUNCT[bp] @ e1) / (1 - PUNCT[bp] @ P)
    w = []
    for y, d in cross_states:
        n = shape_n(y)
        x = (n @ e1) / (1 - n @ P)
        if (x - xa) * (xa - xb) > 0:
            c = "a"
        elif (x - xb) * (xb - xa) > 0:
            c = "b"
        else:
            continue
        w.append(c if d == low else c.upper())
    return "".join(w)


def crossing_states(sol, ts, T):
    out = []
    for t, y in zip(sol.t_events[0], sol.y_events[0]):
        if ts < t <= ts + T:
            out.append((y, 1 if darea(y) > 0 else -1))
    return out
