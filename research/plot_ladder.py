"""Figures for the ladder of periodic orbits A(Ba)^k B a(bA)^k b."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from threebody import Solver  # noqa: E402

BG, PANEL, INK, MUTED, GRID = "#0a0e1a", "#0f1526", "#e6e9f2", "#8a93ab", "#1f2840"
BODY = ["#c2841a", "#3a91d4", "#de5577"]
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": PANEL, "axes.edgecolor": GRID, "axes.labelcolor": MUTED,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "font.size": 10,
                     "axes.titlesize": 11, "axes.titleweight": "bold", "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6})


def ic(m3, v1, v2):
    m = np.array([1.0, 1.0, m3])
    x = np.array([[-1.0, 0, 0], [1.0, 0, 0], [0, 0, 0]])
    v = np.array([[v1, v2, 0], [v1, v2, 0], [-2 * v1 / m3, -2 * v2 / m3, 0]])
    return m, x, v


def orbit(m3, v1, v2, T, n=4000):
    m, x, v = ic(m3, v1, v2)
    return Solver(m, x, v, tol=1e-13).integrate(T, dt_out=T / n)


claims = json.load(open(HERE / "results/claims.json"))
by = {c["id"]: c for c in claims}
rows05 = [by["catalogue I.A i.c. 2 (0.5)"], by["catalogue I.A i.c. 3 (0.5)"], by["L3@0.5"], by["L4@0.5"], by["L5@0.5"]]

fig, axes = plt.subplots(1, 5, figsize=(17, 4.4))
fig.subplots_adjust(left=0.02, right=0.99, top=0.8, bottom=0.04, wspace=0.08)
fig.suptitle("A ladder of periodic three-body orbits, m1 = m2 = 1, m3 = 0.5", x=0.02, ha="left", fontsize=15, fontweight="bold")
fig.text(0.02, 0.86, "word A(Ba)^k B a(bA)^k b · k = 1, 2 are in the Li–Jing–Liao (2018) catalogue · k = 3, 4, 5 are not", color=MUTED, fontsize=10)
for ax, c, k in zip(axes, rows05, [1, 2, 3, 4, 5]):
    tr = orbit(0.5, c["v1"], c["v2"], c["T"])
    for i in range(3):
        ax.plot(tr.x[:, i, 0], tr.x[:, i, 1], color=BODY[i], lw=1.0 if i < 2 else 0.8, alpha=0.95)
        ax.scatter(*tr.x[0, i, :2], s=26 if i < 2 else 14, color=BODY[i], edgecolor=BG, zorder=5)
    ax.set_aspect("equal", adjustable="datalim"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    status = "catalogued" if k <= 2 else "new"
    ax.set_title(f"k = {k}  ({status})", loc="left", color=INK if k > 2 else MUTED)
    ax.text(0.03, 0.03, f"T = {c['T']:.3f}", transform=ax.transAxes, color=MUTED, fontsize=9)
fig.savefig(HERE / "results/ladder_m05.png", dpi=130)

# period vs k at three masses
fig, ax = plt.subplots(figsize=(7.2, 4.4))
fig.subplots_adjust(left=0.11, right=0.97, top=0.86, bottom=0.13)
series = {
    "m3 = 0.5": [(1, by["catalogue I.A i.c. 2 (0.5)"]["T"]), (2, by["catalogue I.A i.c. 3 (0.5)"]["T"])] + [(c["k"], c["T"]) for c in claims if c["id"].endswith("@0.5")],
    "m3 = 0.25": [(c["k"], c["T"]) for c in claims if c["id"].startswith("L") and c["id"].endswith("@0.25")],
    "m3 = 0.1": [(c["k"], c["T"]) for c in claims if c["id"].endswith("@0.1")],
}
cols = ["#3a91d4", "#c2841a", "#de5577"]
for (name, pts), col in zip(series.items(), cols):
    pts = sorted(pts)
    ks, Ts = zip(*pts)
    ax.plot(ks, Ts, "-", color=col, lw=2)
    ax.scatter(ks, Ts, s=36, color=col, edgecolor=BG, zorder=5)
    ax.text(ks[-1] + 0.2, Ts[-1], name, color=INK, va="center", fontsize=9)
ax.axvspan(0.5, 2.5, color="#8a93ab", alpha=0.08)
ax.text(1.5, 95, "catalogued", color=MUTED, ha="center", fontsize=9)
ax.set_xlabel("rung k"); ax.set_ylabel("period T")
ax.set_xlim(0.5, 13); ax.set_ylim(0, 105)
ax.set_title("Each rung adds one repetition and a constant time", loc="left")
fig.savefig(HERE / "results/ladder_periods.png", dpi=130)
print("ok")
