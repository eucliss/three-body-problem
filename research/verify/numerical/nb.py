"""Independent planar Newtonian 3-body integrator (G=1), written from scratch.

State Y = [x1,y1,x2,y2,x3,y3, vx1,vy1,vx2,vy2,vx3,vy3].
"""
import numpy as np
from scipy.integrate import solve_ivp


def make_rhs(m):
    m1, m2, m3 = m

    def rhs(t, Y):
        x1, y1, x2, y2, x3, y3 = Y[0], Y[1], Y[2], Y[3], Y[4], Y[5]
        dx12 = x2 - x1; dy12 = y2 - y1
        dx13 = x3 - x1; dy13 = y3 - y1
        dx23 = x3 - x2; dy23 = y3 - y2
        r12 = (dx12 * dx12 + dy12 * dy12) ** -1.5
        r13 = (dx13 * dx13 + dy13 * dy13) ** -1.5
        r23 = (dx23 * dx23 + dy23 * dy23) ** -1.5
        ax1 = m2 * dx12 * r12 + m3 * dx13 * r13
        ay1 = m2 * dy12 * r12 + m3 * dy13 * r13
        ax2 = -m1 * dx12 * r12 + m3 * dx23 * r23
        ay2 = -m1 * dy12 * r12 + m3 * dy23 * r23
        ax3 = -m1 * dx13 * r13 - m2 * dx23 * r23
        ay3 = -m1 * dy13 * r13 - m2 * dy23 * r23
        return np.array([Y[6], Y[7], Y[8], Y[9], Y[10], Y[11],
                         ax1, ay1, ax2, ay2, ax3, ay3])
    return rhs


def ic(m3, v1, v2):
    Y = np.zeros(12)
    Y[0:2] = (-1.0, 0.0)
    Y[2:4] = (1.0, 0.0)
    Y[4:6] = (0.0, 0.0)
    Y[6:8] = (v1, v2)
    Y[8:10] = (v1, v2)
    Y[10:12] = (-2.0 * v1 / m3, -2.0 * v2 / m3)
    return Y


def energy(Y, m):
    p = Y[:6].reshape(3, 2); v = Y[6:].reshape(3, 2)
    K = 0.5 * sum(m[i] * v[i] @ v[i] for i in range(3))
    U = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            U -= m[i] * m[j] / np.linalg.norm(p[i] - p[j])
    return K + U


def angmom(Y, m):
    p = Y[:6].reshape(3, 2); v = Y[6:].reshape(3, 2)
    return sum(m[i] * (p[i, 0] * v[i, 1] - p[i, 1] * v[i, 0]) for i in range(3))


def momentum(Y, m):
    v = Y[6:].reshape(3, 2)
    return sum(m[i] * v[i] for i in range(3))


def pair_dists(P):
    """P: (..., 6) positions -> (..., 3) distances r12, r13, r23"""
    x1, y1, x2, y2, x3, y3 = [P[..., k] for k in range(6)]
    return np.stack([np.hypot(x2 - x1, y2 - y1), np.hypot(x3 - x1, y3 - y1),
                     np.hypot(x3 - x2, y3 - y2)], axis=-1)


def integrate(m, Y0, T, rtol=1e-13, atol=1e-13, dense=True, method='DOP853'):
    sol = solve_ivp(make_rhs(m), (0.0, T), Y0, method=method, rtol=rtol,
                    atol=atol, dense_output=dense)
    return sol
