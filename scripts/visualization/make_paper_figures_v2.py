#!/usr/bin/env python3
"""
Paper figures, rebuilt on the 2026-09-26 canonical sources (best-NLL run of ten).
2026-09-30: sources switched to the accurate-integral canonical fits (tau <= 0 set to 0)
and the stage-6 split-half reruns; Fig. 3 is a single panel (PBCM removed from the paper).

WHY A NEW FILE.  make_paper_figures.py reads sources that are all superseded:
the shipped-density fits (refit_gateB_*), the break-bug 4c file
(fit_4c_static_nDT3_*), the 4i decomposition (dropped from the model line on
2026-08-11), and the pre-correction split-half run.  It also draws panels the
manuscript no longer reports: a non-decision-time axis (a ceiling artefact) and
a 4i comparison forest.

F3  model comparison        onset is close to fit-neutral; the contrast weight
                            alone buys nothing; the bootstrap null sits where
                            the observed difference sits
F4  group differences       three axes: common-signal drift, threshold, onset
F5  split-half              direction agreement only, no effect sizes
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RUN  = ROOT/"results"/"tDDM"
OUT  = ROOT/"docs"/"PAPER"/"figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "axes.linewidth": 0.8,
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    # Springer Nature asks for Helvetica or Arial.  Arial draws U+03C4 as a bare
    # T that reads as a Latin letter, and matplotlib registers only the regular
    # face of Helvetica.ttc, so bold would silently render unbold under a
    # Helvetica-first stack.  Text is therefore Arial and tau alone comes from
    # Helvetica through mathtext.  Same convention as make_schematic_figures.py.
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Helvetica",
    "mathtext.it": "Helvetica", "pdf.fonttype": 42, "ps.fonttype": 42,
})
HC_C, DEP_C = "#4C72B0", "#C44E52"
NULL_C, OBS_C = "#BBBBBB", "#333333"

# ---- canonical sources (2026-08-20) -----------------------------------------
f4j = pd.read_csv(RUN/"canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv")
f4j["tau"] = f4j["tau"].clip(lower=0)          # every tau <= 0 is the contrast model
f4a = pd.read_csv(RUN/"corrected_likelihood_20260819_012638/fits_4a_static1.csv")
f4c = pd.read_csv(RUN/"corrected_likelihood_20260819_005659/fits_4c_static.csv")
# mean-first and contrast-first split-half runs, fitted separately under the
# same reduced budget, concatenated so F5 keeps both orderings.
shs  = pd.concat([
    pd.read_csv(RUN/"split_half_predictive_image_20260928_115546/split_half_summary.csv"),
    pd.read_csv(RUN/"split_half_predictive_image_20260928_134940/split_half_summary.csv")],
    ignore_index=True)

def panel_label(ax, letter, dx=-0.16, dy=1.06):
    """Lower-case bold panel label, as Communications Psychology requires for
    multi-part figures ("labelled with a lower-case bold a, b, and so on")."""
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left")


def split_g(df, col):
    return df.loc[df.group == "HC", col].values, df.loc[df.group == "DEP", col].values

def cohen_d(a, b):
    na, nb = len(a), len(b)
    s = np.sqrt(((na-1)*a.var(ddof=1) + (nb-1)*b.var(ddof=1)) / (na+nb-2))
    d = (b.mean() - a.mean()) / s
    se = np.sqrt((na+nb)/(na*nb) + d*d/(2*(na+nb-2)))
    return d, d-1.96*se, d+1.96*se

def fmt(v, sign=True):
    """Format with the Unicode minus so labels match the tick formatter."""
    t = f"{v:+.2f}" if sign else f"{v:.2f}"
    return t.replace("-", "\u2212")

def pstr(p):
    return "p < .001" if p < .001 else f"p = {p:.3f}".replace("0.", ".")

def stars(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "n.s."

# =============================================================================
# F3 — model comparison
# =============================================================================
def fig_model_comparison():
    fig, ax = plt.subplots(figsize=(3.6, 3.5))

    # (a) BIC relative to the baseline model.  Plotting absolute BIC needs a
    # truncated axis (the three values sit within one unit of 75), which
    # exaggerates the differences; differences from the baseline start at zero.
    d_ca = float(np.mean(f4c.BIC.values - f4a.BIC.values))
    d_ja = float(np.mean(f4j.BIC.values - f4a.BIC.values))
    d_jc = float(np.mean(f4j.BIC.values - f4c.BIC.values))
    names = ["baseline\n(k = 3)", "contrast\n(k = 4)", "time-contrast\n(k = 5)"]
    vals = [0.0, d_ca, d_ja]
    bars = ax.bar(names, vals, color=["#DDDDDD", "#BBBBBB", "#8899BB"],
                  edgecolor="k", linewidth=0.7, width=0.62)
    ax.axhline(0, color="k", lw=0.8)
    # 2026-09-30: all three differences are now >= 0, so the axis starts just below zero and the
    # contrast-to-time-contrast step is drawn above the bars (it was drawn below when it was negative)
    span = max(abs(min(vals)), abs(max(vals)))
    ax.set_ylim(-span*0.12, span*1.62)
    for b, v in zip(bars[1:], vals[1:]):
        ax.text(b.get_x()+b.get_width()/2, v+span*0.04, fmt(v), ha="center", va="bottom", fontsize=8.5)
    ya = span*1.30
    ax.annotate("", xy=(1.0, ya), xytext=(2.0, ya), arrowprops=dict(arrowstyle="<->", lw=0.8, color="k"))
    ax.annotate(fmt(d_jc), xy=(1.5, ya + span*0.05), ha="center", va="bottom", fontsize=8.5)
    ax.set_ylabel("mean ΔBIC per participant\n(relative to baseline)")
    ax.set_title("Nested model comparison", loc="left")

    fig.savefig(OUT/"F3_model_comparison.png", bbox_inches="tight")
    fig.savefig(OUT/"F3_model_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  F3  dBIC contrast-baseline {d_ca:+.3f} | time-contrast-baseline {d_ja:+.3f} | time-contrast-contrast {d_jc:+.3f}")

# =============================================================================
# F4 — group differences on three axes
# =============================================================================
def violin_pair(ax, hc, dep, title, ylabel, highlight=False):
    parts = ax.violinplot([hc, dep], positions=[0, 1], widths=0.78,
                          showextrema=False, showmedians=False)
    for pc, c in zip(parts["bodies"], [HC_C, DEP_C]):
        pc.set_facecolor(c); pc.set_alpha(0.30); pc.set_edgecolor(c); pc.set_linewidth(0.9)
    bp = ax.boxplot([hc, dep], positions=[0, 1], widths=0.20, patch_artist=True,
                    showfliers=False, medianprops=dict(color="k", lw=1.1),
                    whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8))
    for patch, c in zip(bp["boxes"], [HC_C, DEP_C]):
        patch.set_facecolor("white"); patch.set_edgecolor(c); patch.set_linewidth(1.0)
    ax.plot([0, 1], [hc.mean(), dep.mean()], "D", ms=4.5, color="k", zorder=5)
    p = stats.ttest_ind(dep, hc, equal_var=False).pvalue
    d, _, _ = cohen_d(hc, dep)
    ylo = min(hc.min(), dep.min()); yhi = max(hc.max(), dep.max())
    rng = yhi - ylo
    ax.set_ylim(ylo - rng*0.08, yhi + rng*0.30)
    ybar = yhi + rng*0.11
    ax.plot([0, 0, 1, 1], [ybar, ybar+rng*0.035, ybar+rng*0.035, ybar], lw=0.9, color="k")
    ax.text(0.5, ybar + rng*0.055, stars(p), ha="center", va="bottom", fontsize=9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["HC", "DEP"])
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold" if highlight else "normal")
    if highlight:
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_linewidth(1.3)

def fig_group_differences():
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.2))
    specs = [("d_mean", "Mean-signal drift", "drift rate", False),
             ("thres",  "Threshold", "threshold", False),
             ("tau",    "Estimated onset", "$\\tau$ (s)", True)]
    for ax, (col, title, yl, hl), lab in zip(axes, specs, "abc"):
        hc, dep = split_g(f4j, col)
        violin_pair(ax, hc, dep, title, yl, highlight=hl)
        panel_label(ax, lab)
    fig.subplots_adjust(wspace=0.42)
    fig.savefig(OUT/"F4_group_differences.png", bbox_inches="tight")
    fig.savefig(OUT/"F4_group_differences.pdf", bbox_inches="tight")
    plt.close(fig)
    for col in ("d_mean", "thres", "tau"):
        hc, dep = split_g(f4j, col)
        print(f"  F4  {col:8s} d={cohen_d(hc,dep)[0]:+.3f} p={stats.ttest_ind(dep,hc,equal_var=False).pvalue:.4f}")

# =============================================================================
# F5 — split-half direction agreement
# =============================================================================
def fig_split_half():
    order = ["4a", "4c_mean_static", "4j_mean_first", "4j_contrast_first"]
    labels = {"4a": "baseline\ndrift", "4c_mean_static": "contrast model\ndrift",
              "4j_mean_first": "time-contrast $\\tau$\n(mean first)",
              "4j_contrast_first": "time-contrast $\\tau$\n(contrast first)"}
    d = shs.set_index("model")
    order = [m for m in order if m in d.index]
    vals = [d.loc[m, "reliability_agree"]*100 for m in order]
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    bars = ax.bar(range(len(order)), vals, color="#8899BB", edgecolor="k",
                  linewidth=0.7, width=0.6)
    K = int(shs.K.iloc[0])
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, v+1.4, f"{v:.1f}%\n({round(v*K/100)}/{K})",
                ha="center", va="bottom", fontsize=8)
    ax.axhline(50, color="k", lw=0.8, ls=(0, (3, 3)))
    ax.set_xlim(-0.55, len(order) - 1 + 0.85)
    ax.text(len(order) - 1 + 0.80, 52, "chance", fontsize=7.5, ha="right",
            va="bottom", color="#555555")
    ax.set_ylim(0, 122)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([labels[m] for m in order])
    ax.set_ylabel("splits agreeing in direction (%)")
    ax.set_title("Direction agreement across 30 random participant splits", loc="left")
    fig.savefig(OUT/"F5_split_half.png", bbox_inches="tight")
    fig.savefig(OUT/"F5_split_half.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  F5  " + " | ".join(f"{m}={v:.1f}%" for m, v in zip(order, vals)))

if __name__ == "__main__":
    print("building figures on the accurate-integral canonical sources (2026-09-30)")
    fig_model_comparison(); fig_group_differences(); fig_split_half()
    print(f"written to {OUT}")
