"""Render the validation figures in figures/ (dark 'night sky' style)."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from threebody import Solver, analysis, leapfrog, physics, presets, rk4  # noqa: E402

BG, PANEL, INK, MUTED, GRID = "#0a0e1a", "#0f1526", "#e6e9f2", "#8a93ab", "#1f2840"
BODY = ["#c2841a", "#3a91d4", "#de5577"]   # validated categorical steps (dark surface)
GLOW = ["#f2b84b", "#7cc4f5", "#f58ea6"]   # brighter tints for the star cores

plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": PANEL, "axes.edgecolor": GRID,
    "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
})


def fading_trail(ax, xy, color, lw=1.6):
    pts = xy[:, None, :]
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    alpha = np.linspace(0.08, 1.0, len(segs))
    rgba = np.array([matplotlib.colors.to_rgba(color, a) for a in alpha])
    ax.add_collection(LineCollection(segs, colors=rgba, linewidths=lw, capstyle="round"))


def star(ax, p, i, m):
    s = 40 + 40 * m ** (2 / 3)
    ax.scatter(*p, s=s * 6, color=GLOW[i], alpha=0.12, lw=0, zorder=5)
    ax.scatter(*p, s=s, color=GLOW[i], edgecolor=BG, lw=1.5, zorder=6)


def orbit_panel(ax, m, x, v, T, title, sub, frac=1.0):
    tr = Solver(m, x, v, tol=1e-13).integrate(T, dt_out=T / 1500)
    start = int(len(tr.t) * (1 - frac))
    for i in range(3):
        fading_trail(ax, tr.x[start:, i, :2], BODY[i])
        star(ax, tr.x[-1, i, :2], i, m[i])
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale_view()
    ax.margins(0.1)
    ax.grid(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, loc="left")
    ax.text(0.02, 0.03, sub, transform=ax.transAxes, color=MUTED, fontsize=8.5)
    return tr


def main():
    orbits = json.loads((ROOT / "data" / "periodic_orbits.json").read_text())
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.4))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.9, bottom=0.07, wspace=0.18, hspace=0.28)
    fig.suptitle("The three-body problem, solved numerically", x=0.05, ha="left",
                 fontsize=17, fontweight="bold")
    fig.text(0.05, 0.935, "Regularized Bulirsch–Stoer (LogH leapfrog + Richardson extrapolation), tolerance 1e-13, G = 1",
             color=MUTED, fontsize=10)

    m, x, v = presets.figure_eight()
    tr = orbit_panel(axes[0, 0], m, x, v, 6.32591398, "Figure-eight choreography",
                     "Moore 1993 · Chenciner–Montgomery 2000 · trails: last T/3", frac=1/3)
    o = orbits["butterfly_I"]
    m, x, v = presets.suvakov(o["p1"], o["p2"])
    tr = orbit_panel(axes[0, 1], m, x, v, o["T"], "Butterfly I",
                     f"Šuvakov–Dmitrašinović 2013 · Newton-refined, closes to {o['residual']:.0e}")
    o = orbits["moth_I"]
    m, x, v = presets.suvakov(o["p1"], o["p2"])
    tr = orbit_panel(axes[0, 2], m, x, v, o["T"], "Moth I",
                     f"T = {o['T']:.6f} · closes to {o['residual']:.0e}")

    # Pythagorean problem
    ax = axes[1, 0]
    m, x, v = presets.pythagorean()
    tr = Solver(m, x, v, tol=1e-13).integrate(62, dt_out=0.01)
    for i in range(3):
        fading_trail(ax, tr.x[:, i, :2], BODY[i], lw=1.0)
        star(ax, tr.x[-1, i, :2], i, m[i])
        xy = np.clip(tr.x[-1, i, :2], [-5.6, -8.6], [5.6, 6.6])
        ax.annotate(f"m = {int(m[i])}", xy, xytext=(10, -4 + 10 * (i - 1)), textcoords="offset points",
                    color=INK, fontsize=9)
    c = analysis.classify(tr.x[-1], tr.v[-1], m)
    ax.set_aspect("equal", adjustable="datalim"); ax.set_xlim(-6, 6); ax.set_ylim(-9, 7)
    ax.set_title("Burrau's Pythagorean problem (3, 4, 5)", loc="left")
    ax.text(0.02, 0.97, f"t = 62: mass 3 ejected (v∞ = {c['v_inf']:.2f})\nmasses 4 + 5 bound, e = {c['e']:.3f}",
            transform=ax.transAxes, color=MUTED, fontsize=8.5, va="top")

    # energy error vs time for three integrators
    ax = axes[1, 1]
    E0 = physics.energy(x, v, m)
    s = Solver(m, x, v, tol=1e-12)
    ts, es = [], []
    while s.t < 20:
        s.step()
        ts.append(s.t); es.append(abs(physics.energy(s.x, s.v, m) / E0 - 1))
    dt, n = 1e-3, 20000
    for name, fn, col in (("Leapfrog, dt = 1e-3", leapfrog, "#6f7893"), ("RK4, dt = 1e-3", rk4, "#a9b1c7")):
        X, V = x.copy(), v.copy()
        tt, ee = [], []
        for k in range(200):
            _, X, V = fn(m, X, V, dt, n // 200)
            tt.append((k + 1) * dt * n / 200); ee.append(abs(physics.energy(X, V, m) / E0 - 1))
        ax.semilogy(tt, ee, color=col, lw=1.6)
    ax.text(2.3, 3e-2, "RK4 and leapfrog (dt = 1e-3)\nblow up at the first close\nencounter, t ≈ 1.9", color="#a9b1c7", fontsize=8.5)
    ax.semilogy(ts, np.maximum(es, 1e-16), color=INK, lw=2)
    ax.text(19.6, 4e-14, "regularized solver, tol 1e-12", color=INK, ha="right", fontsize=8.5)
    ax.set_xlim(0, 20); ax.set_ylim(1e-16, 10)
    ax.set_xlabel("time"); ax.set_ylabel("|ΔE / E₀|")
    ax.set_title("Energy error through close encounters", loc="left")

    # butterfly effect: two universes 1e-4 apart
    ax = axes[1, 2]
    a = Solver(m, x, v, tol=1e-13).integrate(60, dt_out=0.01)
    x2 = x.copy(); x2[0, 0] += 1e-4
    b = Solver(m, x2, v, tol=1e-13).integrate(60, dt_out=0.01)
    k = np.searchsorted(a.t, 20)
    for i in range(3):
        ax.plot(a.x[k:, i, 0], a.x[k:, i, 1], color=BODY[i], lw=1.2)
        ax.plot(b.x[k:, i, 0], b.x[k:, i, 1], color=BODY[i], lw=1.0, ls=(0, (2, 2)), alpha=0.8)
        star(ax, a.x[-1, i, :2], i, m[i]); star(ax, b.x[-1, i, :2], i, m[i])
    ca = analysis.classify(a.x[-1], a.v[-1], m); cb = analysis.classify(b.x[-1], b.v[-1], m)
    ax.set_aspect("equal", adjustable="datalim"); ax.set_xlim(-7, 7); ax.set_ylim(-9, 7)
    ax.text(0.02, 0.97, f"solid: original → mass {int(m[ca['third']])} escapes\n"
            f"dashed: body 1 moved 10⁻⁴ → mass {int(m[cb['third']])} escapes",
            transform=ax.transAxes, color=MUTED, fontsize=8.5, va="top")
    ax.set_title("Butterfly effect: same problem, nudged 10⁻⁴", loc="left")

    out = ROOT / "figures" / "overview.png"
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main()
