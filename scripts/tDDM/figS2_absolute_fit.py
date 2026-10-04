#!/usr/bin/env python
"""
Supplementary Figure S2 -- the response-time signature of the mean and
contrast axes, with the 4j prediction overlaid.  It is the picture of
Table S9 and of identifiability analysis 5 in Table S6.

Two changes from figures/labfeedback_2026_08/fig3_rt_mean_contrast.png:

  1. WITHIN-SUBJECT centring.  The old panels binned pooled trials, so
     between-participant RT differences (HC ~3.0 s vs DEP ~3.5 s) leaked into
     every curve.  That leak is not cosmetic: with the pooled binning the
     |contrast| slope looks positive for the model too, and only a
     within-subject analysis shows that the model has it backwards.
  2. The MODEL PREDICTION is drawn.  Each participant's own trials are simulated
     from their own fitted 4j parameters (5 draws averaged per trial), so the
     comparison is like for like.

Result: 4j reproduces the signed mean effect but reverses the |contrast| effect
and hugely overstates the |mean| effect -- it produces the textbook diffusion
pattern (slowest when drift is near zero) that the data do not show.

Source: results/tDDM/rt_signature_check_20260928_005349/
Run: conda activate hkusrp && python scripts/tDDM/figS2_absolute_fit.py
"""
import glob
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "results", "tDDM", "rt_signature_check_20260928_005349")
# trial_level.csv holds observed RTs and is not shipped in the public repository;
# use a run directory given on the command line, else the newest run that has it.
if len(sys.argv) > 1:
    SRC = sys.argv[1]
elif not os.path.exists(os.path.join(SRC, "trial_level.csv")):
    runs = sorted(glob.glob(os.path.join(ROOT, "results", "tDDM", "rt_signature_check_*", "trial_level.csv")))
    if not runs:
        sys.exit("run scripts/tDDM/rt_signature_model_check.R first (trial_level.csv not found)")
    SRC = os.path.dirname(runs[-1])
FIG = os.path.join(ROOT, "docs", "PAPER", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200,
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
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
OBS, MOD = "#1b1b1b", "#1e7a56"

d = pd.read_csv(os.path.join(SRC, "trial_level.csv"))
d["abs_mean"] = d.mean_sig.abs()
d["abs_contrast"] = d.contrast_sig.abs()

# within-subject centring:each participant's own mean log RT removed, so the curves
# show deviation from that participant's typical speed
for c in ("rt_obs", "rt_sim"):
    d["l_" + c] = np.log(d[c])
    d["c_" + c] = d["l_" + c] - d.groupby("sub")["l_" + c].transform("mean")


def curve(col, nbin=5):
    """Quintile-bin WITHIN each participant, then average the centred log RT."""
    q = d.groupby("sub")[col].transform(
        lambda v: pd.qcut(v.rank(method="first"), nbin, labels=False))
    g = d.assign(q=q).groupby("q")
    return (g[col].mean(),
            g["c_rt_obs"].mean(), g["c_rt_obs"].sem(),
            g["c_rt_sim"].mean(), g["c_rt_sim"].sem())


PANELS = [("mean_sig", "a", "Common value", "mean signal"),
          ("contrast_sig", "b", "Contrast (signed)", "contrast signal"),
          ("abs_contrast", "c", "Contrast magnitude", "|contrast|")]

fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.3), constrained_layout=True)

for k, (col, tag, title, xlab) in enumerate(PANELS):
    a = ax[k]
    x, mo, so, mm, sm = curve(col)
    a.axhline(0, color="0.75", lw=0.7, zorder=0)
    a.errorbar(x, mo, yerr=1.96 * so, fmt="o-", color=OBS, lw=1.5, ms=4.5,
               capsize=2.5, label="observed")
    a.errorbar(x, mm, yerr=1.96 * sm, fmt="s--", color=MOD, lw=1.5, ms=4,
               capsize=2.5, label="4j prediction")
    a.set_xlabel(xlab)
    if k == 0:
        a.set_ylabel("log RT, within-subject")
    a.text(-0.22, 1.05, tag, transform=a.transAxes, fontsize=11,
           fontweight="bold", va="bottom", ha="left")
    a.set_title(title, loc="left")

# panel d: the mixed-model coefficients, data against model
cmp = pd.read_csv(os.path.join(SRC, "comparison.csv"))
sim = pd.read_csv(os.path.join(SRC, "sim_coefficients.csv"))
TERMS = ["mean_sig", "abs_mean", "contrast_sig", "abs_contrast"]
LAB = ["mean", "|mean|", "contrast", "|contrast|"]
y = np.arange(len(TERMS))[::-1]

a = ax[3]
a.axvline(0, color="0.55", lw=0.8)
obs_b = [cmp.loc[cmp.term == t, "observed"].iloc[0] for t in TERMS]
mod_b = [sim[t].mean() for t in TERMS]
mod_e = [1.96 * sim[t].std(ddof=1) for t in TERMS]
a.errorbar(mod_b, y + .13, xerr=mod_e, fmt="s", color=MOD, ms=4.5, capsize=2.5,
           label="4j prediction")
a.plot(obs_b, y - .13, "o", color=OBS, ms=5, label="observed")
a.set_yticks(y); a.set_yticklabels(LAB)
a.set_xlabel("β on log RT")
a.text(-0.22, 1.05, "d", transform=a.transAxes, fontsize=11,
       fontweight="bold", va="bottom", ha="left")
a.set_title("Mixed model", loc="left")
a.set_ylim(-0.6, len(TERMS) - 0.4)

h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="outside upper center", ncol=2, handlelength=1.8,
           columnspacing=2.0)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"figS2_absolute_fit.{ext}"),
                bbox_inches="tight", pad_inches=0.03)
plt.close(fig)

print("observed vs model (mixed-model β on log RT)")
for t, o, m, e in zip(LAB, obs_b, mod_b, mod_e):
    # a sign difference only counts when the model's interval actually excludes
    # the observed value -- for `contrast` both estimates sit on top of zero
    outside = abs(o - m) > e
    flag = ("  <- SIGN FLIP" if outside and np.sign(o) != np.sign(m)
            else "  <- outside" if outside else "")
    print(f"  {t:12s} obs {o:+.4f}   model {m:+.4f} ± {e:.4f}{flag}")
print("\nsaved:", FIG)
