#!/usr/bin/env python
"""
Supplementary Figure S1 -- recovery of the individual onset.

WHY THIS FILE EXISTS.  The figure was in docs/PAPER/figures/ with no generator
anywhere in the repository, so it could not be reproduced or refonted.  This
script rebuilds it from the canonical CSVs.

2026-08-21.  Two changes on the PI's instruction.
  * The density-validation figure this file also produced was dropped.  Its
    three panels replotted Table S2 and its y axis ran to 140% to fit the
    misspecified density, which compressed the 0.34-5.87% range under test.
  * The misspecified-density series was dropped.  That density was our own
    extension and was never published, so its behaviour is not the reader's
    concern.
  * The panels became two scatters, at 20 and at 80 trials.  The line plots they
    replace put Table S5's r column and Error SD column on axes, which the table
    already carries; the scatters show the shrinkage and the spread instead.

SOURCE
  results/tDDM/recovery_redesign_tau_20260928_192519  (arm redesign)
  r = -0.077 at 20 trials and +0.406 at 80; error SD 1.155 and 0.765.

Run: conda activate hkusrp && python scripts/tDDM/figS1_individual_recovery.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R    = os.path.join(ROOT, "results", "tDDM")
FIG  = os.path.join(ROOT, "docs", "PAPER", "figures")

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 300,
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "legend.frameon": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    # Arial for the lettering, tau alone from Helvetica through mathtext:
    # Arial draws U+03C4 as a bare T.  See docs/PAPER/figures/README.md.
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Helvetica",
    "mathtext.it": "Helvetica", "pdf.fonttype": 42, "ps.fonttype": 42,
})
CORR = "#1e7a56"
TAU = r"$\tau$"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.join(FIG, name))


def panel_title(ax, letter, title, dx=-0.24):
    """Bold letter outside the axes and a plain left title, as F3-F7 use."""
    ax.text(dx, 1.05, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="left")
    ax.set_title(title, loc="left")


# =============================================================================
# S1 -- individual onset recovery under the corrected likelihood
# =============================================================================
def arm(run, what):
    p = os.path.join(R, run, f"recovery_{what}.csv")
    d = pd.read_csv(p)
    return d[d["arm"] == "redesign"]


def figS1():
    summ = arm("recovery_redesign_tau_20260928_192519", "summary").set_index("n_trials")
    fits = arm("recovery_redesign_tau_20260928_192519", "fits")

    fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.9), sharex=True, sharey=True)
    lo = min(fits.true_tau.min(), fits.est_tau.min())
    hi = max(fits.true_tau.max(), fits.est_tau.max())
    pad = 0.06 * (hi - lo)
    grid = np.linspace(lo, hi, 50)

    for ax, n, letter in zip(axes, (20, 80), "ab"):
        f = fits[fits["n_trials"] == n]
        ax.plot(grid, grid, ":", lw=1.1, color="black", label="identity", zorder=1)
        ax.scatter(f.true_tau, f.est_tau, s=11, color=CORR, alpha=0.35,
                   linewidths=0, zorder=2)
        b1, b0 = np.polyfit(f.true_tau, f.est_tau, 1)
        ax.plot(grid, b0 + b1 * grid, lw=1.9, color=CORR, zorder=3,
                label="least squares")
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"generating {TAU} (s)")
        ax.text(0.04, 0.94, f"r = {summ.loc[n, 'r_tau']:+.2f}".replace("-", "\u2212"),
                transform=ax.transAxes,
                va="top", ha="left")
        panel_title(ax, letter, f"{n} trials", dx=-0.20)
    axes[0].set_ylabel(f"recovered {TAU} (s)")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.06))
    fig.subplots_adjust(wspace=0.12)
    save(fig, "figS1_individual_recovery")
    print(f"  S1  r {summ.loc[20,'r_tau']:+.3f} / {summ.loc[80,'r_tau']:+.3f}"
          f"   errSD {summ.loc[20,'sd_err_tau']:.3f} / {summ.loc[80,'sd_err_tau']:.3f}"
          f"   n per panel {len(fits[fits.n_trials==20])}")


if __name__ == "__main__":
    figS1()
