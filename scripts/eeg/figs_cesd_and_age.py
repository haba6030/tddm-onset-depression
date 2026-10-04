#!/usr/bin/env python
"""Figures for the two PI-requested robustness analyses.

figS3_cesd_continuous.png  — CES-D as a continuous predictor: per-parameter scatter
                          with OLS fit, r and p, and the cutoff-16 line that the
                          dichotomised analysis uses.
figS4_age_control.png      — age by cohort, the correlation of age with each fitted
                          parameter, and the group contrast before and after
                          adjustment for age and gender.

Both write to docs/PAPER/figures/. Run from anywhere:
    python scripts/tDDM/figs_cesd_and_age.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# 2026-08-14: switched to the CORRECTED-transition-density fits. The previous
# default (refit_gateB_20260624_172551/fits_4j_widenDT.csv) came from the
# misspecified split-at-tau density; see 260808_LABMEETING_FEEDBACK.md §3-4.
# 2026-09-26: canonical fits are the highest-likelihood of ten optimizer runs.
# DEoptim is stochastic and 4j does not converge in 200 generations, so a single
# run is not reproducible. See docs/PAPER/TODO1_Rerun_Tracker.md.
# 2026-09-30: accurate-integral canonical fits; tau <= 0 set to 0 in load().
FITS = os.path.join(ROOT, "results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv")
ENDTOEND = os.path.join(ROOT, "results/tDDM/end_to_end_4j_nDT3_image_20260810_043550/per_rep_results.csv")
OUTDIR = os.path.join(ROOT, "docs/PAPER/figures")
RT_LO, RT_HI, CUTOFF = 0.15, 10.0, 16

HC_C, DEP_C, LINE_C, GREY = "#3B7DB5", "#D1603D", "#1F1F1F", "#595959"
PARAMS = [("d_mean", "Mean-signal drift"), ("thres", "Threshold"),
          ("nDT", "Non-decision time (s)"), ("tau", "Fitted onset $\\tau$ (s)"),
          ("d_contrast", "Contrast weight")]

plt.rcParams.update({
    # Arial for the lettering, tau alone from Helvetica through mathtext:
    # Arial draws U+03C4 as a bare T.  See make_schematic_figures.py.
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Helvetica",
    "mathtext.it": "Helvetica", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.edgecolor": GREY, "xtick.color": GREY, "ytick.color": GREY,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 11,
})


def panel_letter(ax, letter, dx=-0.14, dy=1.05):
    """Bold lower-case letter outside the axes, the convention of F3-F7."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="left")


def load():
    fits = pd.read_csv(FITS)
    fits["tau"] = fits["tau"].clip(lower=0)          # every tau <= 0 is the contrast model
    frames = []
    for name, grp in [("study4_HC.csv", "HC"), ("study4_dep.csv", "DEP")]:
        d = pd.read_csv(os.path.join(ROOT, "data", name))
        d["grp"] = grp
        frames.append(d)
    allt = pd.concat(frames, ignore_index=True)
    ok = allt["choice_RT"].between(RT_LO, RT_HI, inclusive="right")
    agg = (allt[ok].groupby("sub")
           .agg(grp=("grp", "first"), ces=("cesScore", "first"),
                age=("age", "first"), gender=("gender", "first"),
                medRT=("choice_RT", "median"), n_used=("choice_RT", "size"))
           .reset_index())
    m = fits.merge(agg, left_on="subject", right_on="sub")
    m["dep"] = (m["grp"] == "DEP").astype(float)
    return m


def pstr(p):
    return "p < .001" if p < .001 else f"p = {p:.3f}".replace("0.", ".")


def ols(X, y):
    """Return (betas, ses, ts, ps) for design X (with intercept column)."""
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = b / se
    return b, se, t, 2 * stats.t.sf(np.abs(t), dof)


# ---------------------------------------------------------------- figure 1
def fig_cesd(m):
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.2))
    for letter, (ax, (key, label)) in zip("abcde", zip(axes.ravel(), PARAMS)):
        panel_letter(ax, letter)
        y, x = m[key].values, m["ces"].values.astype(float)
        r, p = stats.pearsonr(x, y)
        for g, c in (("HC", HC_C), ("DEP", DEP_C)):
            s = m["grp"] == g
            ax.scatter(x[s], y[s], s=15, color=c, alpha=.45, edgecolor="none",
                       label=f"{g} (n = {s.sum()})", zorder=2)
        xs = np.linspace(x.min(), x.max(), 60)
        b, se, _, _ = ols(np.column_stack([np.ones_like(x), x]), y)
        ax.plot(xs, b[0] + b[1] * xs, color=LINE_C, lw=2, zorder=4)
        # 95% band
        n = len(x)
        sxx = ((x - x.mean()) ** 2).sum()
        resid = y - (b[0] + b[1] * x)
        s = np.sqrt(resid @ resid / (n - 2))
        half = 1.96 * s * np.sqrt(1 / n + (xs - x.mean()) ** 2 / sxx)
        ax.fill_between(xs, b[0] + b[1] * xs - half, b[0] + b[1] * xs + half,
                        color=LINE_C, alpha=.13, lw=0, zorder=3)
        ax.axvline(CUTOFF, color=GREY, ls=":", lw=1.3, zorder=1)
        ax.set_title(f"{label}\nr = {r:+.3f},  {pstr(p)}",
                     fontweight="bold" if p < .05 else "normal", loc="left",
                     color="black" if p < .05 else GREY)
        ax.set_xlabel("CES-D total")
        ax.set_ylabel(label)
    ax = axes.ravel()[-1]
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    h, l = axes[0, 0].get_legend_handles_labels()
    ax.legend(h, l, fontsize=10, frameon=False, loc="upper left",
              bbox_to_anchor=(0.0, 0.80))
    fig.tight_layout(pad=1.3)
    out = os.path.join(OUTDIR, "figS3_cesd_continuous.png")
    fig.savefig(out, dpi=200, facecolor="white")
    fig.savefig(out.replace(".png", ".pdf"), facecolor="white")
    print("wrote", out)


# ---------------------------------------------------------------- figure 2
SHORT = {"Mean-signal drift": "Drift", "Threshold": "Threshold",
         "Non-decision time (s)": "nDT", "Fitted onset $\\tau$ (s)": "Onset $\\tau$",
         "Contrast weight": "Contrast"}


def fig_age(m):
    fig = plt.figure(figsize=(13.0, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.5],
                          left=.055, right=.985, top=.80, bottom=.145, wspace=.45)

    # (A) age distribution by group
    axA = fig.add_subplot(gs[0])
    bins = np.arange(15, 80, 4)
    for g, c in (("HC", HC_C), ("DEP", DEP_C)):
        a = m.loc[m["grp"] == g, "age"]
        axA.hist(a, bins=bins, alpha=.55, color=c, label=f"{g}: {a.mean():.1f} ({a.std():.1f})")
        axA.axvline(a.mean(), color=c, lw=2)
    t, p = stats.ttest_ind(m.loc[m.grp == "DEP", "age"], m.loc[m.grp == "HC", "age"],
                           equal_var=False)
    d = (m.loc[m.grp == "DEP", "age"].mean() - m.loc[m.grp == "HC", "age"].mean()) / m["age"].std()
    axA.set_xlabel("Age (years)")
    axA.set_ylabel("Participants")
    axA.set_ylim(0, axA.get_ylim()[1] * 1.32)
    axA.legend(fontsize=9.5, frameon=False, loc="upper right")
    panel_letter(axA, "a")

    # (B) why it does not threaten the result: r(age, param)
    axB = fig.add_subplot(gs[1])
    rs = [stats.pearsonr(m["age"].astype(float), m[k])[0] for k, _ in PARAMS]
    ps = [stats.pearsonr(m["age"].astype(float), m[k])[1] for k, _ in PARAMS]
    ypos = np.arange(len(PARAMS))[::-1]
    axB.barh(ypos, rs, color=["#8A4B14" if q < .05 else "#C9C9C9" for q in ps], height=.6)
    axB.axvline(0, color="black", lw=1)
    axB.set_yticks(ypos)
    axB.set_yticklabels([SHORT[lb] for _, lb in PARAMS], fontsize=9.5)
    axB.tick_params(axis="y", pad=2)
    for yy, rr, pp in zip(ypos, rs, ps):
        inside = rr > .10
        axB.text(rr - .008 if inside else rr + .010, yy, f"{rr:+.3f}",
                 va="center", ha="right" if inside else "left", fontsize=8.5,
                 color="white" if inside else "black",
                 fontweight="bold" if pp < .05 else "normal")
    axB.set_xlim(-.09, .30)
    axB.set_xticks([0, .1, .2, .3])
    axB.set_xlabel("r (age, parameter)")
    panel_letter(axB, "b")

    # (C) standardized group effect, unadjusted vs age+gender adjusted
    axC = fig.add_subplot(gs[2])
    rows = []
    for k, lb in PARAMS:
        y = ((m[k] - m[k].mean()) / m[k].std()).values
        n = len(y)
        X1 = np.column_stack([np.ones(n), m["dep"].values])
        X2 = np.column_stack([np.ones(n), m["dep"].values,
                              m["age"].values.astype(float), m["gender"].values.astype(float)])
        b1, se1, _, p1 = ols(X1, y)
        b2, se2, _, p2 = ols(X2, y)
        rows.append((lb, b1[1], se1[1], p1[1], b2[1], se2[1], p2[1]))
    ypos = np.arange(len(rows))[::-1].astype(float)
    for i, (lb, bu, su, pu, ba, sa, pa) in enumerate(rows):
        yy = ypos[i]
        axC.errorbar(bu, yy + .16, xerr=1.96 * su, fmt="o", ms=6, color=GREY,
                     capsize=3, lw=1.5, zorder=3)
        axC.errorbar(ba, yy - .16, xerr=1.96 * sa, fmt="o", ms=6, color="#8A4B14",
                     capsize=3, lw=1.8, zorder=3)
        axC.text(.80, yy + .16, pstr(pu), fontsize=8.5, color=GREY, va="center")
        axC.text(.80, yy - .16, pstr(pa), fontsize=8.5, color="#8A4B14",
                 va="center", fontweight="bold" if pa < .05 else "normal")
    axC.axvline(0, color="black", lw=1)
    axC.set_yticks(ypos)
    axC.set_yticklabels([SHORT[r[0]] for r in rows], fontsize=9.5)
    axC.tick_params(axis="y", pad=2)
    axC.set_ylim(-0.9, len(rows) - 0.4)
    axC.set_xlim(-.72, 1.20)
    axC.set_xticks([-.5, -.25, 0, .25, .5])
    axC.set_xlabel("Standardized group effect (DEP − HC), 95% CI")
    axC.plot([], [], "o", color=GREY, label="unadjusted")
    axC.plot([], [], "o", color="#8A4B14", label="+ age + gender")
    axC.legend(fontsize=9, frameon=False, loc="lower left", ncol=2,
               bbox_to_anchor=(-0.02, -0.03))
    panel_letter(axC, "c")
    out = os.path.join(OUTDIR, "figS4_age_control.png")
    fig.savefig(out, dpi=200, facecolor="white")
    fig.savefig(out.replace(".png", ".pdf"), facecolor="white")
    print("wrote", out)
    for lb, bu, su, pu, ba, sa, pa in rows:
        print(f"  {lb:24s} unadj {bu:+.3f} ({pstr(pu)})   adj {ba:+.3f} ({pstr(pa)})")




# ---------------------------------------------------------------- figure 3
def fig_endtoend():
    """Group-level end-to-end recovery: a CALIBRATION check, not an existence test."""
    r = pd.read_csv(ENDTOEND)
    truth = -r["true_delta"].iloc[0]          # flip to "DEP later = positive"
    rec = -r["recovered_delta"].values
    lo, hi = -r["ci_hi"].values, -r["ci_lo"].values
    mean_rec = rec.mean()

    fig, ax = plt.subplots(figsize=(10.4, 4.0))
    y = np.arange(len(rec))[::-1]
    ax.axvline(0, color="black", lw=1.4)
    ax.axvline(truth, color=DEP_C, ls="--", lw=1.8)
    ax.text(truth + .012, len(rec) - .35,
            f"TRUE Δ$\\tau$ put into the simulation = {truth:.3f} s",
            color=DEP_C, fontsize=9.5, va="center")
    txt_x = max(hi.max(), truth) + .07
    covers = [(l <= truth <= h) for l, h in zip(lo, hi)]
    for i, (yy, m, l, h) in enumerate(zip(y, rec, lo, hi)):
        excl = l > 0 or h < 0
        ax.errorbar(m, yy, xerr=[[m - l], [h - m]], fmt="o", ms=7,
                    color="#8A4B14", capsize=4, lw=1.8, zorder=3)
        tag = ("excludes 0" if excl else "includes 0") + \
              (" · covers truth" if covers[i] else " · MISSES truth")
        ax.text(txt_x, yy, tag, fontsize=9, va="center",
                color="#8A4B14" if excl else GREY,
                fontweight="bold" if excl else "normal")
    ax.axvline(mean_rec, color="#8A4B14", lw=1.2, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([f"replicate {i + 1}" for i in range(len(rec))], fontsize=10)
    ax.set_ylim(-0.7, len(rec) - 0.15)
    ax.set_xlim(min(lo.min(), 0) - .15, txt_x + .52)
    ax.set_xlabel("Recovered Δ$\\hat{\\tau}$  (DEP − HC, s),  95% CI")
    n_sign = int(sum(np.sign(rec) == np.sign(truth)))
    n_excl = int(sum((l > 0) or (h < 0) for l, h in zip(lo, hi)))
    ax.set_title(
        f"Feed a gap in as truth → replicate spread is wide and NOT one-sided\n"
        f"mean recovered {mean_rec:.3f} s = {100 * mean_rec / truth:.0f}% of the injected "
        f"{truth:.3f} s · sign {n_sign}/{len(rec)} · covers truth {sum(covers)}/{len(rec)} · "
        f"excludes 0 {n_excl}/{len(rec)}\n"
        f"CAVEAT: the injected {truth:.3f} s is the PRE-FIX observed gap; the corrected-density gap "
        f"is +0.288 s, so this control sits at the wrong operating point",
        fontsize=9.8, fontweight="bold", loc="left", color="#1F1F1F")
    fig.tight_layout(pad=1.1)
    out = os.path.join(OUTDIR, "archive", "FS_endtoend_group_recovery.png")
    fig.savefig(out, dpi=200, facecolor="white")
    fig.savefig(out.replace(".png", ".pdf"), facecolor="white")
    print("wrote", out)
    print(f"  truth={truth:.4f} recovered={rec.round(4)} mean={mean_rec:.4f} "
          f"({100 * mean_rec / truth:.1f}%) excl0={sum((l > 0) or (h < 0) for l, h in zip(lo, hi))}/{len(rec)}")


if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)
    m = load()
    print(f"N = {len(m)}  (HC {(m.grp == 'HC').sum()}, DEP {(m.grp == 'DEP').sum()})")
    fig_cesd(m)
    fig_age(m)
    # 2026-08-21: the end-to-end control is not one of the six identifiability
    # analyses the Supplement reports (Table S6), and it was run at the pre-fix
    # observed gap of 0.169 s rather than the canonical 0.288 s, which its own
    # caption states.  It is kept for the record and written to archive/ only.
    if "--with-endtoend" in sys.argv:
        fig_endtoend()
