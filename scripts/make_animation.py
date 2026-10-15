"""Render figures/pythagorean.mp4: Burrau's problem from rest to ejection (needs ffmpeg)."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.collections import LineCollection
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from threebody import Solver, analysis, physics, presets  # noqa: E402

BG, INK, MUTED = "#060912", "#e6e9f2", "#949db5"
BODY = ["#c2841a", "#3a91d4", "#de5577"]
GLOW = ["#f6c66a", "#8fd0ff", "#ff9bb3"]


def main(t_end=64.0, fps=30, seconds=24):
    m, x, v = presets.pythagorean()
    n = fps * seconds
    tr = Solver(m, x, v, tol=1e-13).integrate(t_end, dt_out=t_end / n)
    E0 = physics.energy(x, v, m)
    lo = tr.x[:, :, :2].min((0, 1)) - 1.2
    hi = tr.x[:, :, :2].max((0, 1)) + 1.2
    hi[1] += 1.6  # headroom for the title
    # pad the narrower side so the frame is a 4:5 portrait (1080 x 1350)
    c, span = (lo + hi) / 2, hi - lo
    if span[1] / span[0] > 1.25: span[0] = span[1] / 1.25
    else: span[1] = span[0] * 1.25
    lo, hi = c - span / 2, c + span / 2
    fig = plt.figure(figsize=(7.2, 9.0), dpi=150, facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1], facecolor=BG)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal"); ax.axis("off")
    rng = np.random.default_rng(3)
    ax.scatter(rng.uniform(lo[0], hi[0], 320), rng.uniform(lo[1], hi[1], 320), s=rng.uniform(0.2, 2.2, 320) ** 2,
               color="#cfd6ee", alpha=0.35, lw=0)
    trails = [LineCollection([], linewidths=1.6, capstyle="round") for _ in range(3)]
    for lc in trails:
        ax.add_collection(lc)
    halos = [ax.scatter([], [], s=(40 + 30 * mi) * 7, color=GLOW[i], alpha=0.14, lw=0, zorder=4) for i, mi in enumerate(m)]
    cores = [ax.scatter([], [], s=40 + 30 * mi, color=GLOW[i], edgecolor=BG, lw=1.2, zorder=5) for i, mi in enumerate(m)]
    title = ax.text(0.04, 0.965, "Burrau's Pythagorean problem", transform=ax.transAxes, color=INK, fontsize=15, fontweight="bold")
    sub = ax.text(0.04, 0.94, "masses 3, 4, 5 released from rest · G = 1", transform=ax.transAxes, color=MUTED, fontsize=9.5)
    clock = ax.text(0.04, 0.03, "", transform=ax.transAxes, color=INK, fontsize=10, family="monospace")
    status = ax.text(0.04, 0.055, "", transform=ax.transAxes, color=MUTED, fontsize=9.5)
    for i in range(3):
        ax.text(0.70 + 0.09 * i, 0.03, f"m={int(m[i])}", transform=ax.transAxes, color=GLOW[i], fontsize=9.5, family="monospace")
    trail_len = 3.5
    writer = FFMpegWriter(fps=fps, bitrate=6000, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    out = ROOT / "figures" / "pythagorean.mp4"
    with writer.saving(fig, str(out), dpi=150):
        for k in range(len(tr.t)):
            t = tr.t[k]
            lo = np.searchsorted(tr.t, t - trail_len)
            for i in range(3):
                pts = tr.x[lo : k + 1, i, :2]
                if len(pts) > 1:
                    segs = np.stack([pts[:-1], pts[1:]], 1)
                    a = np.linspace(0.0, 1.0, len(segs)) ** 1.4
                    trails[i].set_segments(segs)
                    trails[i].set_color([matplotlib.colors.to_rgba(BODY[i], q) for q in a])
                halos[i].set_offsets([tr.x[k, i, :2]]); cores[i].set_offsets([tr.x[k, i, :2]])
            dE = abs(physics.energy(tr.x[k], tr.v[k], m) / E0 - 1)
            clock.set_text(f"t = {t:6.2f}   |ΔE/E₀| = {max(dE, 1e-16):.1e}")
            c = analysis.classify(tr.x[k], tr.v[k], m)
            status.set_text(f"mass {int(m[c['third']])} escaping · masses {int(m[c['pair'][0]])}+{int(m[c['pair'][1]])} bound, e = {c['e']:.3f}"
                            if c["state"] == "escape" else "three-way interplay: close encounters, no pattern")
            writer.grab_frame()
    print(out)


if __name__ == "__main__":
    main()
