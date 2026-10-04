#!/usr/bin/env python3
"""
F1 and F2 -- the two conceptual schematics of the manuscript.

WHY A SCRIPT AND NOT AN IMAGE MODEL.  Springer Nature does not accept
generative-AI images as published artwork, and the two figures carry structural
constraints an image model does not hold: the phase-2 slope in F1 must equal the
contrast-model slope in F2b, and the three F2 panels must be identical in
geometry so that only the traces differ.

WHAT IS FITTED AND WHAT IS ILLUSTRATIVE.  The geometry follows the canonical
mean-of-10 fits (results/tDDM/canonical_4j_meanof10_20260820_164804):
boundary separation 2*thres with thres = 1.99 (drawn as +/-a with a = 1 in
figure units), onset tau = 0.73 s (canonical mean 0.727).  Three quantities are
drawn for legibility rather than to scale, and the captions say the traces are
illustrative:
  * drift.  The phase-2 slope is drawn at three times the phase-1 slope, not at
    the fitted ratio.  A trial-level drift depends on that trial's signal
    values, and 42.4% of d_contrast estimates sit at the bound, so no single
    fitted pair is representative; the ratio is set so that the change of slope
    at tau, which is what the figure exists to show, survives print reduction.
  * diffusion noise.  Drawn at 0.13, not the engine's sigma = 1.4.  At the true
    sigma a single path swamps the boundary separation, and the ordering the
    three F2 panels must convey inverts from seed to seed.
  * non-decision time.  Drawn as a short segment.  The fitted nDT (~2.1 s) sits
    against a censored ceiling (15.8% at the 3.0 s bound) and drawing it to
    scale would consume the panel.

SEEDS.  Every path is drawn from a fixed seed.  The seeds were screened so that
each panel's path stays close to its own mean drift (peak deviation ~0.1a): at a seed with a large
excursion the F2 ordering inverts on the page and panel b reads as shallower
than panel a, which is the opposite of what the panel is there to show.

IN-PANEL TEXT IS DELIBERATELY MINIMAL.  Definitions, parameter counts and model
names live in the captions in docs/PAPER/v2/02_Methods.md.
F1 uses: Time, Evidence, Accept, Reject, nDT, tau, Phase 1, Phase 2.
F2 uses: a, b, c, +a, -a, tau.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT  = ROOT/"docs"/"PAPER"/"figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 600,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    # Springer Nature asks for Helvetica or Arial.  Neither mathtext route
    # works here: stixsans falls back to STIXGeneral (a serif) for the italic
    # boundary label and encodes tau in a private-use slot, so the PDF extracts
    # as \ue1eb and cannot be searched or reflowed.  Everything is therefore
    # plain Unicode text.  Tau comes from Helvetica because Arial draws U+03C4
    # as a bare T that reads as a Latin letter (checked against Helvetica,
    # Arial Unicode MS, Verdana, Tahoma, Lucida Grande, DejaVu Sans, Geneva).
    "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.8, "lines.solid_capstyle": "round",
})

INK, GREY = "#1A1A1A", "#9A9A9A"
TAU_FONT  = FontProperties(family="Helvetica", size=8)
BOUND     = 1.0          # boundary in figure units (+a / -a)
TAU       = 0.73         # s, canonical mean onset
NDT       = 0.38         # s, illustrative (see docstring)
V1, V2    = 0.26, 0.95   # drift in figure units per second (phase 2 = 3.7x phase 1)
NOISE     = 0.13

def path(v_pre, v_post, tau, t_end, seed, dt=0.002):
    """One diffusion path; drift switches from v_pre to v_post at tau."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, t_end + dt, dt)
    v = np.where(t < tau, v_pre, v_post)
    x = np.concatenate([[0.0], np.cumsum(v[:-1]*dt + NOISE*np.sqrt(dt)*rng.standard_normal(len(t)-1))])
    hit = np.argmax(x >= BOUND) if (x >= BOUND).any() else None
    if hit:
        t, x = t[:hit+1], x[:hit+1]
        x[-1] = BOUND
    return t, x

def frame(ax, xmax, tick_labels=True):
    """Shared panel furniture: two boundaries, start dot, time axis at zero."""
    for y in (BOUND, -BOUND):
        ax.plot([-NDT, xmax], [y, y], color=INK, lw=1.1, zorder=2)
    ax.plot([-NDT, xmax], [0, 0], color=GREY, lw=0.7, zorder=1)
    ax.plot([-NDT, 0], [0, 0], color=INK, lw=1.4, zorder=3)
    ax.plot([0], [0], "o", color=INK, ms=3.2, zorder=4)
    ax.set_xlim(-NDT-0.16, xmax)
    ax.set_ylim(-BOUND-0.62, BOUND+0.62)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

def bracket(ax, x0, x1, y, label, down=False, tip=0.075, pad=0.055):
    """Square bracket spanning x0..x1 with a centred label."""
    s = -1 if down else 1
    ax.plot([x0, x0, x1, x1], [y, y+s*tip, y+s*tip, y], color=INK, lw=0.7,
            clip_on=False, zorder=5)
    ax.text((x0+x1)/2, y+s*(tip+pad), label, ha="center",
            va="bottom" if not down else "top", color=INK, clip_on=False)

# ---------------------------------------------------------------- F1 ---------
def figure1():
    fig, ax = plt.subplots(figsize=(3.46, 2.20))
    t, x = path(V1, V2, TAU, 3.0, seed=9)
    t_hit, xmax = t[-1], t[-1] + 0.30
    frame(ax, xmax)

    ax.plot([TAU, TAU], [-BOUND-0.44, BOUND+0.10], ls=(0, (3, 2.4)), lw=0.8, color=INK, zorder=2)
    ax.plot(t, x, color=INK, lw=1.0, zorder=4)

    ax.text(-NDT+0.02, BOUND+0.07, "Accept", ha="left", va="bottom", color=INK)
    ax.text(-NDT+0.02, -BOUND-0.07, "Reject", ha="left", va="top", color=INK)
    ax.text(TAU, -BOUND-0.52, "τ", ha="center", va="top", color=INK, fontproperties=TAU_FONT)
    bracket(ax, -NDT, 0, -0.30, "nDT", down=True)
    bracket(ax, 0, TAU, BOUND+0.14, "Phase 1")
    bracket(ax, TAU, t_hit, BOUND+0.14, "Phase 2")
    ax.plot([t_hit], [BOUND], "o", color=INK, ms=3.2, zorder=5)

    ax.annotate("", xy=(xmax, -BOUND-0.44), xytext=(-NDT-0.16, -BOUND-0.44),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK,
                                mutation_scale=7, shrinkA=0, shrinkB=0))
    ax.annotate("", xy=(-NDT-0.16, BOUND+0.30), xytext=(-NDT-0.16, -BOUND-0.44),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK,
                                mutation_scale=7, shrinkA=0, shrinkB=0))
    ax.text(xmax, -BOUND-0.36, "Time", ha="right", va="bottom", color=INK)
    ax.text(-NDT-0.28, 0, "Evidence", rotation=90, ha="right", va="center", color=INK)

    fig.subplots_adjust(left=0.13, right=0.985, top=0.90, bottom=0.13)
    for ext in ("png", "pdf"):
        fig.savefig(OUT/f"F1_model_schematic.{ext}")
    return fig

# ---------------------------------------------------------------- F2 ---------
def figure2():
    xmax = 2.05
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 1.85))
    traces = [path(V1, V1, xmax, xmax, seed=6),
              path(V2, V2, xmax, xmax, seed=19),
              path(V1, V2, TAU, xmax, seed=9)]
    for k, (ax, (t, x)) in enumerate(zip(axes, traces)):
        frame(ax, xmax)
        if k == 2:
            ax.plot([TAU, TAU], [-BOUND, BOUND], ls=(0, (3, 2.4)), lw=0.8, color=INK, zorder=2)
            ax.text(TAU, -BOUND-0.12, "τ", ha="center", va="top", color=INK, fontproperties=TAU_FONT)
        if k == 1:
            ax.plot([0, xmax], [0, V1*xmax], ls=(0, (2.6, 2.2)), lw=0.8, color=GREY, zorder=3)
        ax.plot(t, x, color=INK, lw=1.0, zorder=4)
        if x[-1] >= BOUND:
            ax.plot([t[-1]], [BOUND], "o", color=INK, ms=3.2, zorder=5)
        ax.text(-NDT+0.02, BOUND+0.07, "+a", ha="left", va="bottom", color=INK, style="italic")
        ax.text(-NDT+0.02, -BOUND-0.07, "−a", ha="left", va="top", color=INK, style="italic")
        ax.text(-NDT-0.16, BOUND+0.42, "abc"[k], ha="left", va="bottom",
                fontweight="bold", fontsize=9)

    fig.subplots_adjust(left=0.015, right=0.995, top=0.84, bottom=0.08, wspace=0.10)
    for ext in ("png", "pdf"):
        fig.savefig(OUT/f"F2_model_set.{ext}")
    return fig

if __name__ == "__main__":
    plt.close(figure1()); plt.close(figure2())
    print("wrote", OUT/"F1_model_schematic.png", OUT/"F2_model_set.png")
