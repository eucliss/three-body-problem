"""Validation of the three-body solver against exact solutions and published results."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from threebody import Solver, analysis, leapfrog, physics, presets, rk4  # noqa: E402


def run(preset, t, tol=1e-12, **kw):
    m, x, v = preset
    s = Solver(m, x, v, tol=tol, **kw)
    return m, x, v, s.integrate(t)


def test_ten_integrals_conserved_through_close_encounters():
    m, x, v, tr = run(presets.pythagorean(), 70.0)
    assert tr.stats["dE_rel"] < 1e-9
    I0 = physics.integrals(x, v, m)
    I1 = physics.integrals(tr.x[-1], tr.v[-1], m)
    assert np.allclose(I1["P"], I0["P"], atol=1e-11)
    assert np.allclose(I1["L"], I0["L"], atol=1e-9)
    assert np.allclose(I1["com"], I0["com"], atol=1e-9)


def test_pythagorean_problem_outcome():
    # Szebehely & Peters (1967): mass 3 is ejected, masses 4 and 5 form an eccentric binary.
    m, x, v, tr = run(presets.pythagorean(), 80.0)
    c = analysis.classify(tr.x[-1], tr.v[-1], m)
    assert c["state"] == "escape"
    assert m[c["third"]] == 3
    assert sorted(m[list(c["pair"])]) == [4, 5]
    assert 0.95 < c["e"] < 1.0


def test_matches_scipy_dop853():
    m, x, v = presets.pythagorean()
    t_end = 15.0
    tr = Solver(m, x, v, tol=1e-13).integrate(t_end)

    def f(t, y):
        X = y[:9].reshape(3, 3)
        return np.concatenate([y[9:], physics.accelerations(X, m).ravel()])

    ref = solve_ivp(f, (0, t_end), np.concatenate([x.ravel(), v.ravel()]),
                    method="DOP853", rtol=1e-13, atol=1e-13)
    # chaotic system with close encounters: agreement to ~1e-6 means both are right
    assert np.abs(tr.x[-1].ravel() - ref.y[:9, -1]).max() < 1e-5


def test_figure_eight_is_periodic():
    m, x, v, tr = run(presets.figure_eight(), 6.32591398)
    assert np.abs(tr.x[-1] - x).max() < 1e-6  # limited by the 8-digit published ICs
    assert tr.stats["dE_rel"] < 1e-12


@pytest.mark.parametrize("masses", [(1, 1, 1), (1, 2, 3), (5, 0.1, 0.01)])
def test_lagrange_triangle_rotates_rigidly(masses):
    m, x, v, tr = run(presets.lagrange(masses), 3.0)
    d0 = [np.linalg.norm(x[i] - x[j]) for i, j in ((0, 1), (1, 2), (0, 2))]
    d1 = [np.linalg.norm(tr.x[-1][i] - tr.x[-1][j]) for i, j in ((0, 1), (1, 2), (0, 2))]
    assert np.allclose(d0, 1.0) and np.allclose(d1, 1.0, atol=1e-8)


def test_routh_criterion():
    assert not presets.routh_stable((1, 1, 1))
    assert presets.routh_stable((1, 0.01, 0.0001))


@pytest.mark.parametrize("masses", [(1, 1, 1), (3, 1, 2), (1, 0.001, 0.5)])
def test_euler_collinear_solution(masses):
    (m, x, v), lam = presets.euler(masses)
    # λ is the positive real root of Euler's quintic
    assert abs(np.polyval(presets.euler_quintic(masses), lam)) < 1e-9 * max(masses) * 10
    s = Solver(m, x, v, tol=1e-13)
    tr = s.integrate(1.0)
    r = lambda X: (np.linalg.norm(X[0] - X[1]), np.linalg.norm(X[1] - X[2]))
    assert np.allclose(r(tr.x[-1]), r(x), rtol=1e-7)
    # still collinear
    X = tr.x[-1]
    assert np.linalg.norm(np.cross(X[1] - X[0], X[2] - X[0])) < 1e-7


def test_refined_periodic_orbits_close():
    import json
    path = Path(__file__).resolve().parents[1] / "data" / "periodic_orbits.json"
    orbits = json.loads(path.read_text())
    for name in ("butterfly_I", "moth_I", "figure_eight_SD"):
        o = orbits[name]
        m, x, v = presets.suvakov(o["p1"], o["p2"])
        tr = Solver(m, x, v, tol=1e-13).integrate(o["T"])
        assert np.abs(tr.x[-1] - x).max() < 1e-9, name


def test_regularization_survives_head_on_collision():
    # two bodies falling straight at each other: a true collision singularity at finite t.
    m = np.array([1.0, 1.0, 1e-3])
    x = np.array([[-0.5, 0, 0], [0.5, 0, 0], [0, 5, 0]], float)
    v = np.zeros((3, 3))
    s = Solver(m, x, v, tol=1e-12)
    tr = s.integrate(3.0)  # free-fall time for this pair is pi/4 * sqrt(1/2) ≈ 0.555
    assert np.isfinite(tr.x).all()
    assert tr.stats["dE_rel"] < 1e-8


def test_regularized_beats_fixed_step_methods():
    m, x, v = presets.pythagorean()
    E0 = physics.energy(x, v, m)
    _, xl, vl = leapfrog(m, x, v, 1e-3, 3000)
    _, xr, vr = rk4(m, x, v, 1e-3, 3000)
    tr = Solver(m, x, v).integrate(3.0)
    err = lambda X, V: abs(physics.energy(X, V, m) / E0 - 1)
    assert tr.stats["dE_rel"] < err(xl, vl) and tr.stats["dE_rel"] < err(xr, vr)


def test_shape_sphere_landmarks():
    m, x, _ = presets.lagrange((1, 1, 1))
    assert abs(abs(physics.shape_sphere(x, m)[2]) - 1) < 1e-12  # pole
    (m, x, _), _ = presets.euler((1, 1, 1))
    assert abs(physics.shape_sphere(x, m)[2]) < 1e-12  # equator


def test_lyapunov_positive_for_chaos():
    m, x, v = presets.pythagorean()
    a = Solver(m, x, v).integrate(10, dt_out=0.1)
    x2 = x.copy()
    x2[0, 0] += 1e-10
    b = Solver(m, x2, v).integrate(10, dt_out=0.1)
    d = [analysis.phase_distance(a.x[i], a.v[i], b.x[i], b.v[i]) for i in range(len(a.t))]
    assert analysis.lyapunov_fit(a.t, d) > 0.1
