"""Initial conditions: exact solutions, periodic choreographies and astrophysical triples.

All in units with G = 1 unless stated. Every function returns (m, x, v) with
x, v of shape (3, 3), already shifted to the barycentric frame.
"""
from __future__ import annotations

import numpy as np

from .physics import to_com_frame


def _pack(m, x, v):
    m = np.asarray(m, float)
    x = np.array(x, float)
    v = np.array(v, float)
    if x.shape[1] == 2:
        x = np.hstack([x, np.zeros((len(m), 1))])
        v = np.hstack([v, np.zeros((len(m), 1))])
    x, v = to_com_frame(x, v, m)
    return m, x, v


# --- exact solutions -----------------------------------------------------------

def figure_eight():
    """Chenciner–Montgomery (2000) figure-eight choreography, equal masses. T ≈ 6.32591398."""
    x1 = np.array([-0.97000436, 0.24308753])
    v3 = np.array([-0.93240737, -0.86473146])
    return _pack([1, 1, 1], [x1, -x1, [0, 0]], [-v3 / 2, -v3 / 2, v3])


def lagrange(m=(1.0, 1.0, 1.0), a=1.0, G=1.0, e=0.0):
    """Lagrange's equilateral solution (1772) for ANY masses: the bodies sit on an
    equilateral triangle that rotates rigidly. e > 0 gives the homographic (pulsating)
    version: every body follows a Kepler ellipse of eccentricity e about the barycentre.
    Linearly stable only if 27(m1 m2 + m2 m3 + m3 m1) < (m1 + m2 + m3)² (Routh 1875)."""
    m = np.asarray(m, float)
    ang = np.radians([90, 210, 330])
    x = a / np.sqrt(3) * np.stack([np.cos(ang), np.sin(ang), np.zeros(3)], 1)
    xc = (m[:, None] * x).sum(0) / m.sum()
    x -= xc
    # start at pericentre of the homographic Kepler motion: speed scales by sqrt(1+e)
    w = np.sqrt(G * m.sum() / a**3) * np.sqrt(1 + e)
    v = w * np.stack([-x[:, 1], x[:, 0], np.zeros(3)], 1)
    return _pack(m, x, v)


def routh_stable(m):
    m1, m2, m3 = m
    return 27 * (m1 * m2 + m2 * m3 + m3 * m1) < (m1 + m2 + m3) ** 2


def euler(m=(1.0, 1.0, 1.0), a=2.0, G=1.0):
    """Euler's collinear solution (1767) for ANY masses, bodies ordered 1-2-3 on a line
    that rotates rigidly. The spacing ratio λ = |x3 - x2| / |x2 - x1| is the single real
    positive root of Euler's quintic; here it is found by bisection on the equivalent
    rigid-rotation condition. `a` is the total length |x3 - x1|."""
    m1, m2, m3 = map(float, m)

    def acc(lam):
        a1 = m2 + m3 / (1 + lam) ** 2
        a2 = -m1 + m3 / lam**2
        a3 = -m1 / (1 + lam) ** 2 - m2 / lam**2
        return a1, a2, a3

    def f(lam):
        a1, a2, a3 = acc(lam)
        return (a2 - a1) - (a3 - a2) / lam

    lo, hi = 1e-6, 1e6
    flo = f(lo)
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        fm = f(mid)
        if np.sign(fm) == np.sign(flo):
            lo, flo = mid, fm
        else:
            hi = mid
    lam = np.sqrt(lo * hi)
    a1, a2, _ = acc(lam)
    omega2 = -(a2 - a1) * G  # spacing x2 - x1 = 1 in these units
    xs = np.array([0.0, 1.0, 1.0 + lam])
    s = a / (1.0 + lam)
    xs *= s
    omega = np.sqrt(omega2 / s**3)
    mm = np.array([m1, m2, m3])
    xs -= (mm * xs).sum() / mm.sum()
    x = np.stack([xs, np.zeros(3), np.zeros(3)], 1)
    v = np.stack([np.zeros(3), omega * xs, np.zeros(3)], 1)
    return _pack(mm, x, v), lam


def euler_quintic(m):
    """Coefficients of Euler's quintic in λ (highest power first), for checking."""
    m1, m2, m3 = m
    return [m1 + m2, 3 * m1 + 2 * m2, 3 * m1 + m2, -(m2 + 3 * m3), -(2 * m2 + 3 * m3), -(m2 + m3)]


# Šuvakov & Dmitrašinović (2013, PRL 110, 114301) periodic orbits:
# x1 = (-1, 0), x2 = (1, 0), x3 = (0, 0); v1 = v2 = (p1, p2), v3 = -2 (p1, p2); m = 1, G = 1.
SUVAKOV = {
    "butterfly_I": (0.30689, 0.12551, 6.2356),
    "butterfly_II": (0.39295, 0.09758, 7.0039),
    "bumblebee": (0.18428, 0.58719, 63.5345),
    "moth_I": (0.46444, 0.39606, 14.8939),
    "moth_II": (0.43917, 0.45297, 28.6703),
    "butterfly_III": (0.40592, 0.23016, 13.8658),
    "moth_III": (0.38344, 0.37736, 25.8406),
    "goggles": (0.08330, 0.12789, 10.4668),
    "butterfly_IV": (0.350112, 0.07934, 79.4759),
    "dragonfly": (0.08058, 0.58884, 21.2710),
    "yarn": (0.55906, 0.34919, 55.5018),
    "yin_yang_Ia": (0.51394, 0.30474, 17.3284),
    "yin_yang_Ib": (0.28270, 0.32721, 10.9626),
    "figure_eight_SD": (0.347111, 0.532728, 6.3259),
}


def suvakov(p1, p2):
    p = np.array([p1, p2])
    return _pack([1, 1, 1], [[-1, 0], [1, 0], [0, 0]], [p, p, -2 * p])


# --- astrophysical and classic chaotic problems --------------------------------

def pythagorean():
    """Burrau's problem (1913): masses 3, 4, 5 at rest on the vertices of a 3-4-5 right
    triangle, each opposite the side of its own length. Szebehely & Peters (1967) showed
    it ends with the lightest body ejected and the other two in a tight binary."""
    return _pack([3, 4, 5], [[1, 3], [-2, -1], [1, -1]], np.zeros((3, 2)))


def kepler_pair(m1, m2, a, e, G=1.0, phase="peri"):
    """Positions/velocities of a two-body orbit about its barycentre, in the xy plane."""
    M = m1 + m2
    if phase == "peri":
        r = a * (1 - e)
        vrel = np.sqrt(G * M * (1 + e) / (a * (1 - e)))
    else:
        r = a * (1 + e)
        vrel = np.sqrt(G * M * (1 - e) / (a * (1 + e)))
    x1 = np.array([-m2 / M * r, 0, 0])
    x2 = np.array([m1 / M * r, 0, 0])
    v1 = np.array([0, -m2 / M * vrel, 0])
    v2 = np.array([0, m1 / M * vrel, 0])
    return x1, v1, x2, v2


def kepler16():
    """Kepler-16 (AB)b, the first known circumbinary planet ("Tatooine").
    Units: solar mass, AU, and time unit 1/(2π) yr so that G = 1.
    Doyle et al. (2011): M_A = 0.6897, M_B = 0.20255, a_AB = 0.22431 AU, e_AB = 0.15944;
    planet 0.333 M_Jup = 3.18e-4 M_sun on a = 0.7048 AU, e ≈ 0.0069."""
    mA, mB, mp = 0.6897, 0.20255, 3.18e-4
    xA, vA, xB, vB = kepler_pair(mA, mB, 0.22431, 0.15944)
    Mb = mA + mB
    ap = 0.7048
    vp = np.sqrt((Mb + mp) / ap)
    xp = np.array([0, ap, 0])
    vpv = np.array([-vp, 0, 0])
    return _pack([mA, mB, mp], [xA, xB, xp], [vA, vB, vpv])


def sitnikov(e=0.2, z0=1.0, vz0=0.0, m3=1e-5):
    """Sitnikov problem: two equal stars on an eccentric orbit, a third (light) body
    oscillating along the axis through their barycentre. For e > 0 the phase space is mixed:
    small oscillations are regular, and Moser proved chaotic orbits exist."""
    x1, v1, x2, v2 = kepler_pair(0.5, 0.5, 1.0, e)
    return _pack([0.5, 0.5, m3], [x1, x2, [0, 0, z0]], [v1, v2, [0, 0, vz0]])


def kozai(inc_deg=80.0):
    """Hierarchical triple with a highly inclined outer companion: the inner orbit
    trades inclination for eccentricity (von Zeipel–Lidov–Kozai cycles)."""
    m1, m2, m3 = 1.0, 0.001, 1.0
    x1, v1, x2, v2 = kepler_pair(m1, m2, 1.0, 0.05)
    ao = 5.0
    vo = np.sqrt((m1 + m2 + m3) / ao)
    i = np.radians(inc_deg)
    # outer body on a circular orbit about the inner barycentre, tilted by i about x
    x3 = np.array([ao, 0.0, 0.0])
    v3 = np.array([0.0, vo * np.cos(i), vo * np.sin(i)])
    return _pack([m1, m2, m3], [x1, x2, x3], [v1, v2, v3])


def trojan(mu=0.001, offset=0.01):
    """Sun–Jupiter–Trojan: a test asteroid near L4, 60° ahead of the planet."""
    m1, m2, m3 = 1.0 - mu, mu, 1e-9
    w = 1.0
    x1 = np.array([-mu, 0, 0])
    x2 = np.array([1 - mu, 0, 0])
    L4 = np.array([0.5 - mu, np.sqrt(3) / 2, 0]) * (1 + offset)
    v = lambda p: w * np.array([-p[1], p[0], 0])
    return _pack([m1, m2, m3], [x1, x2, L4], [v(x1), v(x2), v(L4)])


def free_fall(px, py):
    """Free-fall (Agekyan–Anosova) problem: equal masses released from rest, bodies 1
    and 2 at (∓0.5, 0) and body 3 at (px, py)."""
    return _pack([1, 1, 1], [[-0.5, 0], [0.5, 0], [px, py]], np.zeros((3, 2)))
