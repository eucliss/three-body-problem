"""Newtonian gravity for N point masses (N = 3 is the three-body problem).

The general equation of motion, for every body i:

    d²r_i/dt² = -G * Σ_{j≠i} m_j (r_i - r_j) / (|r_i - r_j|² + ε²)^{3/2}

with ε = 0 for true point masses (ε > 0 "softens" the 1/r singularity).

Conventions used throughout:
    x : (N, 3) positions
    v : (N, 3) velocities
    m : (N,)   masses
    U : the *positive* potential magnitude  U = G Σ_{i<j} m_i m_j / r_ij
    E : total energy  E = T - U
"""
from __future__ import annotations

import numpy as np


def accel_potential(x, m, G=1.0, eps2=0.0):
    """Accelerations a (N,3) and potential magnitude U for all bodies."""
    d = x[np.newaxis, :, :] - x[:, np.newaxis, :]  # d[i, j] = x_j - x_i
    r2 = np.einsum("ijk,ijk->ij", d, d) + eps2
    np.fill_diagonal(r2, 1.0)
    inv_r = 1.0 / np.sqrt(r2)
    np.fill_diagonal(inv_r, 0.0)
    a = G * np.einsum("ijk,j,ij->ik", d, m, inv_r**3)
    U = 0.5 * G * np.einsum("i,j,ij->", m, m, inv_r)
    return a, U


def accelerations(x, m, G=1.0, eps2=0.0):
    return accel_potential(x, m, G, eps2)[0]


def kinetic(v, m):
    return 0.5 * float(np.sum(m * np.einsum("ij,ij->i", v, v)))


def potential(x, m, G=1.0, eps2=0.0):
    return accel_potential(x, m, G, eps2)[1]


def energy(x, v, m, G=1.0, eps2=0.0):
    return kinetic(v, m) - potential(x, m, G, eps2)


def momentum(v, m):
    return (m[:, None] * v).sum(0)


def angular_momentum(x, v, m):
    return (m[:, None] * np.cross(x, v)).sum(0)


def center_of_mass(x, v, m):
    M = m.sum()
    return (m[:, None] * x).sum(0) / M, (m[:, None] * v).sum(0) / M


def to_com_frame(x, v, m):
    """Shift to the barycentric frame: total momentum zero, centre of mass at the origin."""
    xc, vc = center_of_mass(x, v, m)
    return x - xc, v - vc


def integrals(x, v, m, G=1.0, eps2=0.0):
    """The ten classical first integrals (energy, momentum, angular momentum, CoM motion)."""
    xc, vc = center_of_mass(x, v, m)
    return {
        "E": energy(x, v, m, G, eps2),
        "P": momentum(v, m),
        "L": angular_momentum(x, v, m),
        "com": xc,
        "vcom": vc,
    }


def two_body_elements(xi, vi, mi, xj, vj, mj, G=1.0):
    """Relative two-body orbit of the pair (i, j): semi-major axis a, eccentricity e,
    specific orbital energy eps, inclination (deg, relative to the z axis)."""
    r = xj - xi
    w = vj - vi
    mu = G * (mi + mj)
    rn = np.linalg.norm(r)
    eps = 0.5 * w @ w - mu / rn
    h = np.cross(r, w)
    hn = np.linalg.norm(h)
    e_vec = np.cross(w, h) / mu - r / rn
    e = float(np.linalg.norm(e_vec))
    a = -mu / (2 * eps) if eps < 0 else np.inf
    inc = float(np.degrees(np.arccos(np.clip(h[2] / hn, -1, 1)))) if hn > 0 else 0.0
    return {"a": a, "e": e, "eps": eps, "inc": inc, "r": rn}


def shape_sphere(x, m, normal=None):
    """Project the triangle formed by the three bodies onto Montgomery's shape sphere.

    Uses mass-weighted Jacobi vectors. Returns a unit vector (w1, w2, w3):
      * the equator (w3 = 0) holds every collinear configuration, including the Euler
        solutions and the three binary-collision points,
      * the two Lagrange equilateral triangles sit at w3 = ±sqrt(3 m1 m2 m3 M)/(m1m2+m2m3+m3m1),
        which is a pole (|w3| = 1) only when the masses are equal.
    """
    m1, m2, m3 = m
    rho1 = x[1] - x[0]
    rho2 = x[2] - (m1 * x[0] + m2 * x[1]) / (m1 + m2)
    mu1 = m1 * m2 / (m1 + m2)
    mu2 = (m1 + m2) * m3 / (m1 + m2 + m3)
    z1 = np.sqrt(mu1) * rho1
    z2 = np.sqrt(mu2) * rho2
    if normal is None:
        normal = np.array([0.0, 0.0, 1.0])
    w1 = z1 @ z1 - z2 @ z2
    w2 = 2 * z1 @ z2
    w3 = 2 * np.cross(z1, z2) @ normal
    n = np.sqrt(w1 * w1 + w2 * w2 + w3 * w3)
    return np.array([w1, w2, w3]) / n if n > 0 else np.zeros(3)
