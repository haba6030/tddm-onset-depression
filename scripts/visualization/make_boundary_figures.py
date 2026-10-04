#!/usr/bin/env python
"""
F6 and F7 -- the two boundary figures of the Results.

Manuscript_Restructure_Plan_260818.md section 11 asks for two figures that draw
the limits of the onset claim: (a) an ANCOVA ladder with the null control, and
(b) a sensitivity ladder with the inversion.  The plan places both in R3; the
finished draft moved that content, so F6 supports R4 (what the onset absorbs)
and F7 supports R5 (what the design can resolve).

  F6 a  group onset coefficient under increasing control  -> R4 paragraph 3
     b  recovered gap when the true gap is zero           -> R4 paragraph 2
  F7 a  recovered gap against an injected gap             -> R5 paragraph 4
     b  sensitivity and bias against trials per participant
     c  the injected gap implied by the observed one      -> R5 paragraph 4

SIGN.  Everything is the depressed-range minus low-symptom difference, the
convention of the Results and of F4.  The stored null and positive controls hold
the opposite sign, so they are negated on read.

SOURCES (2026-09-30: accurate-integral canonical fits and stage-6 controls)
  canonical fits  results/tDDM/canonical_4j_accurate_bestof10/  (tau <= 0 set to 0)
  null control    results/tDDM/nulltau_corrected_image_20260928_211818/  (8 reps)
  positive control results/tDDM/positive_control_tau_20260928_231637/ (20 trials, 4 levels)
                  results/tDDM/positive_control_tau_20260929_012719/ (80 trials, 2 levels)
  superseded:     canonical_4j_meanof10_20260820_164804/
  null control    results/tDDM/nulltau_corrected_image_20260923_220453/  (8 reps)
  positive control results/tDDM/positive_control_tau_20260924_110653/
  trials ladder   results/tDDM/sensitivity_ladder_summary_20260926.csv, built by
                  scripts/tDDM/collect_sensitivity_ladder.py from the four
                  positive-control runs.  That script documents how trials per
                  participant were attributed, which the run directories do not
                  record.

Run: conda activate hkusrp && python scripts/visualization/make_boundary_figures.py
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R    = ROOT/"results"/"tDDM"
OUT  = ROOT/"docs"/"PAPER"/"figures"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300,
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "legend.frameon": False,
    # Arial for the lettering, tau alone from Helvetica through mathtext.
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Helvetica",
    "mathtext.it": "Helvetica", "pdf.fonttype": 42, "ps.fonttype": 42,
})
OBS_C, NULL_C, TAU_C, REF = "#2b2b2b", "#C44E52", "#1e7a56", "#8e8e8e"
MED_C = "#1b4f72"          # the inverted estimate, and the level means in F7a
OBS_D, OBS_SE = 0.209, 0.104      # accurate-integral canonical group contrast (tau <= 0 set to 0)
RT_LO, RT_HI = 0.15, 10.0         # analysis window, as in the Methods

fits = pd.read_csv(R/"canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv")
fits["tau"] = fits["tau"].clip(lower=0)


def panel_label(ax, letter, dx=-0.20, dy=1.05):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT/f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("written", OUT/name)


# ---------------------------------------------------------------- ANCOVA -----
def ancova_ladder():
    tr = pd.concat([pd.read_csv(ROOT/"data"/"study4_HC.csv"),
                    pd.read_csv(ROOT/"data"/"study4_dep.csv")])
    sub = "sub" if "sub" in tr.columns else "subject"
    ok = tr.choice_RT.between(RT_LO, RT_HI, inclusive="right")
    med = tr[ok].groupby(sub).choice_RT.median().rename("medRT")
    m = fits.join(med, on="subject")
    m["grp"] = (m.group == "DEP").astype(int)
    steps = [("group only",      "tau ~ grp"),
             ("+ median RT",     "tau ~ grp + medRT"),
             ("+ threshold",     "tau ~ grp + medRT + thres"),
             ("+ drift, nDT",    "tau ~ grp + medRT + thres + d_mean + nDT")]
    rows = []
    for lab, f in steps:
        r = smf.ols(f, m).fit()
        ci = r.conf_int().loc["grp"]
        rows.append((lab, r.params["grp"], r.bse["grp"], ci[0], ci[1],
                     r.pvalues["grp"], int(r.df_resid)))
    return pd.DataFrame(rows, columns=["step", "b", "se", "lo", "hi", "p", "df"])


# ------------------------------------------------------------------- F6 ------
def figF6():
    lad = ancova_ladder()
    nl = pd.read_csv(R/"nulltau_corrected_image_20260928_211818/per_rep.csv")
    null = -nl.recovered_dtau.values                      # -> DEP minus HC
    m, se = null.mean(), null.std(ddof=1)/np.sqrt(len(null))

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0))

    ax = axes[0]
    y = np.arange(len(lad))[::-1]
    ax.axvline(0, color="black", lw=0.9)
    ax.errorbar(lad.b, y, xerr=[lad.b - lad.lo, lad.hi - lad.b], fmt="o",
                ms=6, color=OBS_C, capsize=3, lw=1.4)
    ax.set_yticks(y); ax.set_yticklabels(lad.step)
    ax.set_ylim(-0.6, len(lad) - 0.4)
    ax.set_xlabel("group onset difference (s)")
    ax.set_title("Under increasing control", loc="left")
    panel_label(ax, "a", dx=-0.42)

    ax = axes[1]
    ax.axhline(0, color="black", lw=0.9)
    ax.axhspan(m - 1.96*se, m + 1.96*se, color=NULL_C, alpha=0.15)
    ax.axhline(m, color=NULL_C, lw=1.6, label="no true gap")
    ax.plot(np.arange(1, len(null)+1), null, "o", ms=5, color=NULL_C)
    ax.axhline(OBS_D, color=OBS_C, ls="--", lw=1.3, label="observed")
    ax.set_xticks(range(1, len(null)+1, 2))
    ax.set_ylim(-0.02, 0.42)
    ax.set_xlabel("replication"); ax.set_ylabel("recovered gap (s)")
    ax.legend(loc="upper left", handlelength=1.4)
    ax.set_title("With no true gap", loc="left")
    panel_label(ax, "b", dx=-0.24)

    fig.subplots_adjust(wspace=0.55)
    save(fig, "F6_onset_controls")
    for _, r in lad.iterrows():
        print(f"  F6a  {r.step:14s} b={r.b:+.3f} SE={r.se:.3f} "
              f"t({r.df})={r.b/r.se:.2f} p={r.p:.4f}")
    print(f"  F6b  null mean {m:+.3f} (SD {null.std(ddof=1):.3f}, {len(null)} reps) "
          f"= {100*m/OBS_D:.0f}% of the observed {OBS_D:.3f}")


# ------------------------------------------------------------------- F7 ------
def figF7():
    pc = pd.read_csv(R/"positive_control_tau_20260928_231637/per_point.csv")
    X = np.vstack([np.ones(len(pc)), pc.true_gap.values]).T
    beta = np.linalg.lstsq(X, pc.recovered_gap.values, rcond=None)[0]
    res = pc.recovered_gap.values - X @ beta
    V = (res @ res / (len(res) - 2)) * np.linalg.inv(X.T @ X)

    def design(pp):
        d = pd.read_csv(pp); Xd = np.vstack([np.ones(len(d)), d.true_gap.values]).T
        b = np.linalg.lstsq(Xd, d.recovered_gap.values, rcond=None)[0]; e = d.recovered_gap.values - Xd @ b
        Vd = (e @ e / (len(e) - 2)) * np.linalg.inv(Xd.T @ Xd); se = np.sqrt(np.diag(Vd))
        from scipy import stats as st
        tq = st.t.ppf(.975, len(e) - 2)
        return b[1], b[1] - tq*se[1], b[1] + tq*se[1], b[0], se[0]
    d20 = design(R/"positive_control_tau_20260928_231637/per_point.csv")    # four levels
    d80 = design(R/"positive_control_tau_20260929_012719/per_point.csv")    # two levels
    tr_n = np.array([20, 80])
    sl, slo, shi = np.array([d20[0], d80[0]]), np.array([d20[1], d80[1]]), np.array([d20[2], d80[2]])
    bi, bse = np.array([d20[3], d80[3]]), np.array([d20[4], d80[4]])

    rng = np.random.default_rng(0)
    B = 300_000
    draw = rng.multivariate_normal(beta, V, B)
    obs = rng.normal(OBS_D, OBS_SE, B)
    ok = draw[:, 1] > 0            # same rule as scripts/tDDM/memo_0919/stage6_summary_accurate.R
    tru = (obs[ok] - draw[ok, 0]) / draw[ok, 1]
    q = np.percentile(tru, [2.5, 50, 97.5])

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6))

    ax = axes[0]
    gaps = sorted(pc.true_gap.unique())
    for g in gaps:
        v = pc[np.isclose(pc.true_gap, g)].recovered_gap
        ax.plot(np.full(len(v), g) + rng.uniform(-.018, .018, len(v)), v, "o",
                color="#9fc3d6", ms=3.2, alpha=.85, zorder=1)
        ax.plot([g], [v.mean()], "o", color="#1b4f72", ms=6, zorder=3)
    xs = np.array([0, 0.73])
    ax.plot(xs, beta[0] + beta[1]*xs, color=TAU_C, lw=1.8, zorder=2, label="fitted")
    ax.plot(xs, beta[0] + xs, color=REF, ls=":", lw=1.3, label="one to one")
    ax.set_xticks([0, 0.2, 0.4, 0.6])   # injected values (0, 0.085, 0.169, 0.34) overlapped as tick labels
    ax.set_xlim(-0.05, 0.75); ax.set_ylim(-0.02, 0.92)
    ax.set_xlabel("injected gap (s)"); ax.set_ylabel("recovered gap (s)")
    ax.legend(loc="upper left", handlelength=1.4)
    ax.set_title("Response to a true gap", loc="left")
    panel_label(ax, "a", dx=-0.30)

    ax = axes[1]
    ax.axhline(1.0, color=REF, ls=":", lw=1.1)
    ax.axhline(0, color="black", lw=0.9)
    ax.errorbar(tr_n, sl, yerr=[sl-slo, shi-sl], fmt="o-", color=TAU_C, lw=1.6,
                ms=5, capsize=2.5, label="sensitivity")
    ax.errorbar(tr_n, bi, yerr=1.96*bse, fmt="s--", color=NULL_C, lw=1.4, ms=4.5,
                capsize=2.5, label="bias (s)")
    ax.set_xscale("log"); ax.set_xticks(tr_n); ax.set_xticklabels(tr_n)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(15, 105); ax.set_ylim(-0.04, 1.28)
    ax.set_xlabel("trials / participant"); ax.set_ylabel("sensitivity  /  bias")
    ax.legend(loc="upper left", handlelength=1.4)
    ax.set_title("With more trials", loc="left")
    panel_label(ax, "b", dx=-0.28)

    ax = axes[2]
    YT = 0.70
    ax.hist(tru[(tru > -1.5) & (tru < 3.0)], bins=80, color=TAU_C, alpha=.65)   # display range only
    ax.axvline(0, color="black", lw=0.9, ymax=YT)
    ax.axvline(OBS_D, color=OBS_C, ls="--", lw=1.4, ymax=YT, label="observed")
    ax.axvline(q[1], color=MED_C, lw=1.8, ymax=YT, label="implied")
    for v in (q[0], q[2]):
        ax.axvline(v, color=MED_C, ls=(0, (2, 1.6)), lw=1.1, ymax=YT)
    ax.set_xlim(-1.5, 3.0); ax.set_ylim(0, 1.42*ax.get_ylim()[1]); ax.set_yticks([])
    ax.set_xlabel("implied true gap (s)"); ax.set_ylabel("draws")
    ax.legend(loc="upper right", handlelength=1.3)
    ax.set_title("Inverting the response", loc="left")
    panel_label(ax, "c", dx=-0.16)

    fig.subplots_adjust(wspace=0.42)
    save(fig, "F7_onset_sensitivity")
    print(f"  F7a  slope {beta[1]:+.3f}  intercept {beta[0]:+.3f}")
    print(f"  F7c  implied {q[1]:+.3f} s  95% [{q[0]:+.3f}, {q[2]:+.3f}]  "
          f"P(>0) = {(tru > 0).mean():.3f}")


if __name__ == "__main__":
    figF6(); figF7()
