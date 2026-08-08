"""Create publication-style main and supplementary figures from audited outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PIPELINE = ROOT
OUTPUTS = ROOT / "outputs"
FIGURES = ROOT / "paper" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

COLORS = {
    "QMP": "#0072B2",
    "Row-closed": "#D55E00",
    "CLR": "#009E73",
    "CLR (multiplicative)": "#CC79A7",
    "Prevalence": "#7F7F7F",
    "Clinical": "#6A3D9A",
    "Ink": "#20242A",
    "Muted": "#6B7280",
    "Light": "#EEF2F5",
}


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#4B5563",
            "axes.linewidth": 0.7,
            "xtick.color": "#374151",
            "ytick.color": "#374151",
            "text.color": COLORS["Ink"],
            "savefig.dpi": 300,
            "figure.dpi": 120,
        }
    )


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def figure1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def box(x, y, w, h, title, lines, edge, fill="#FFFFFF", title_size=10.0, body_size=8.0):
        patch = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.015",
            linewidth=1.2, edgecolor=edge, facecolor=fill,
        )
        ax.add_patch(patch)
        ax.text(x + 0.02, y + h - 0.032, title, fontsize=title_size, fontweight="bold", color=edge, va="top")
        ax.text(x + 0.02, y + h - 0.073, "\n".join(lines), fontsize=body_size, va="top", linespacing=1.38)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13,
                                     linewidth=1.1, color="#64748B"))

    ax.text(0.02, 0.975, "A", fontsize=12, fontweight="bold", va="top")
    ax.text(0.5, 0.965, "Paired abundance-representation benchmark", ha="center",
            fontsize=12.5, fontweight="bold")
    box(0.04, 0.73, 0.43, 0.18, "LCPM", [
        "589 public profiles",
        "Primary CRC contrast: n=252",
        "205 controls; 47 CRC",
        "676 species; cells per gram",
    ], COLORS["QMP"], "#F2F8FC", body_size=7.6)
    box(0.53, 0.73, 0.43, 0.18, "MetaCardis", [
        "1,087 unique IDs; 994 profiles",
        "Primary disjoint contrast: n=672",
        "369 MMC; 303 IHD",
        "729 MGS; load-corrected index",
    ], COLORS["Clinical"], "#F7F3FA", body_size=7.6)

    arrow(0.255, 0.73, 0.36, 0.62)
    arrow(0.745, 0.73, 0.64, 0.62)
    box(0.12, 0.45, 0.76, 0.17, "Matched representations", [
        "QMP: published quantitative values (cohort-specific units)",
        "Row-closed: deterministic closure of QMP (not native sequencing RMP)",
        "CLR: minimum-positive primary rule + multiplicative sensitivity",
    ], "#374151", "#F8FAFC", body_size=7.8)

    arrow(0.35, 0.45, 0.23, 0.34)
    arrow(0.50, 0.45, 0.50, 0.34)
    arrow(0.65, 0.45, 0.77, 0.34)
    box(0.02, 0.15, 0.30, 0.19, "Feature inference", [
        "Outcome-blind pooled filters",
        "LCPM: 93; MetaCardis: 404 taxa",
        "Benjamini–Hochberg FDR",
    ], COLORS["Row-closed"], "#FFF7F3", title_size=8.8, body_size=6.8)
    box(0.35, 0.15, 0.30, 0.19, "Disease discrimination", [
        "Repeated 5-fold CV ×10",
        "Identical folds; local filtering",
        "Fixed gradient boosting",
    ], COLORS["QMP"], "#F2F8FC", title_size=8.8, body_size=6.8)
    box(0.68, 0.15, 0.30, 0.19, "Robustness", [
        "Source-aligned filter sensitivity",
        "Core vs medication adjustment",
        "49 shared species; no pooling",
    ], COLORS["CLR"], "#F2FAF7", title_size=8.8, body_size=6.8)
    ax.text(0.5, 0.07,
            "Question: how do representation, filtering, zero handling, and covariates change inference and discrimination?",
            ha="center", fontsize=8.2, fontweight="bold")
    save(fig, "Figure_1_study_design")


def figure2_prediction() -> None:
    repeats = pd.read_csv(OUTPUTS / "prediction/cv_repeat_metrics.csv")
    reps = ["QMP", "Row-closed", "CLR", "CLR (multiplicative)"]
    display = {
        "QMP": "QMP",
        "Row-closed": "Row-closed",
        "CLR": "CLR\nmin-positive",
        "CLR (multiplicative)": "CLR\nmultiplicative",
    }
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), sharey=True)
    rng = np.random.default_rng(531)
    for ax, cohort, letter in zip(axes, ["LCPM", "MetaCardis"], ["A", "B"]):
        subset = repeats.loc[(repeats["cohort"] == cohort) & repeats["model"].isin(reps)]
        pivot = subset.pivot(index="repeat", columns="model", values="roc_auc")[reps]
        for _, row in pivot.iterrows():
            ax.plot(range(4), row.values, color="#CBD5E1", linewidth=0.8, alpha=0.75, zorder=1)
        for pos, rep in enumerate(reps):
            values = pivot[rep].to_numpy()
            jitter = rng.uniform(-0.045, 0.045, len(values))
            ax.scatter(pos + jitter, values, s=20, color=COLORS[rep], alpha=0.75,
                       edgecolor="white", linewidth=0.4, zorder=2)
            ax.errorbar(pos, values.mean(), yerr=values.std(ddof=1), fmt="o", ms=6,
                        color=COLORS[rep], ecolor=COLORS[rep], capsize=3,
                        linewidth=1.4, zorder=3)
            ax.text(pos, 0.735, f"{values.mean():.3f}", ha="center", va="bottom", fontsize=8.0,
                    color=COLORS[rep], fontweight="bold")
        ax.axhline(0.5, color="#9CA3AF", linestyle="--", linewidth=0.8)
        ax.set_xticks(range(4), [display[r] for r in reps], fontsize=7.6)
        ax.set_ylim(0.49, 0.76)
        ax.set_title(f"{cohort} primary contrast")
        ax.set_ylabel("Repeated-CV ROC AUC" if ax is axes[0] else "")
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.6)
        panel_letter(ax, letter)
    fig.suptitle("QMP showed no consistent advantage; CLR performance depended on zero replacement", y=1.02,
                 fontsize=11.5, fontweight="bold")
    fig.tight_layout()
    save(fig, "Figure_2_prediction")


def figure3_effect_concordance() -> None:
    lcpm = pd.read_csv(OUTPUTS / "associations/lcpm_crc_vs_ctl_associations.csv")
    lcpm = lcpm.loc[lcpm["filter_specification"].eq("pooled_outcome_blind")]
    lq = lcpm.loc[lcpm["component"].eq("qmp")].set_index("feature")
    lr = lcpm.loc[lcpm["component"].eq("row_closed")].set_index("feature")
    l = lq[["effect", "significant_q_lt_0_05", "q_value_bh"]].join(
        lr[["effect", "significant_q_lt_0_05", "q_value_bh"]], lsuffix="_qmp", rsuffix="_row"
    )

    meta = pd.read_csv(
        OUTPUTS / "associations/metacardis_ihd_vs_mmc_hurdle_qmp_row_closed_clr.csv"
    )
    core = meta.loc[meta["adjustment"].eq("core")]
    mq = core.loc[core["component"].eq("qmp_nonzero")].set_index("matrix_column")
    mr = core.loc[core["component"].eq("row_closed_nonzero")].set_index("matrix_column")
    m = mq[["species", "effect", "significant_q_lt_0_05", "q_value_bh"]].join(
        mr[["effect", "significant_q_lt_0_05", "q_value_bh"]], rsuffix="_row"
    ).rename(columns={"effect": "effect_qmp", "significant_q_lt_0_05": "significant_qmp",
                      "q_value_bh": "q_qmp", "effect_row": "effect_row",
                      "significant_q_lt_0_05_row": "significant_row", "q_value_bh_row": "q_row"})

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.45))
    panels = [
        (
            axes[0], l["effect_qmp"], l["effect_row"],
            l["significant_q_lt_0_05_qmp"] | l["significant_q_lt_0_05_row"],
            "LCPM: CRC vs CTL (outcome-blind filter)",
        ),
        (
            axes[1], m["effect_qmp"], m["effect_row"],
            m["significant_qmp"] | m["significant_row"],
            "MetaCardis: IHD vs MMC (core-adjusted)",
        ),
    ]
    for ax, x, y, sig, title in panels:
        r = float(x.corr(y))
        direction = f"{100 * (np.sign(x) == np.sign(y)).mean():.1f}%"
        if ax is axes[0]:
            qmp_calls = int(l["significant_q_lt_0_05_qmp"].sum())
            row_calls = int(l["significant_q_lt_0_05_row"].sum())
        else:
            qmp_calls = int(m["significant_qmp"].sum())
            row_calls = int(m["significant_row"].sum())
        calls = f"{qmp_calls} / {row_calls}"
        ax.scatter(x[~sig], y[~sig], s=12, color="#9CA3AF", alpha=0.55, linewidth=0)
        ax.scatter(x[sig], y[sig], s=25, color=COLORS["QMP"], alpha=0.9,
                   edgecolor="white", linewidth=0.4)
        lo = min(float(x.min()), float(y.min()))
        hi = max(float(x.max()), float(y.max()))
        pad = (hi - lo) * 0.06
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#64748B",
                linestyle="--", linewidth=0.9)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.axhline(0, color="#E5E7EB", linewidth=0.7)
        ax.axvline(0, color="#E5E7EB", linewidth=0.7)
        ax.set_xlabel("QMP effect")
        ax.set_ylabel("Row-closed effect")
        ax.set_title(title)
        ax.text(0.04, 0.96, f"Pearson r = {r:.3f}\nDirection agreement = {direction}\nFDR calls QMP / row-closed = {calls}",
                transform=ax.transAxes, va="top", fontsize=8.2,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#D1D5DB", alpha=0.9))
    panel_letter(axes[0], "A")
    panel_letter(axes[1], "B")
    fig.suptitle("QMP and row-closed effect patterns were highly concordant",
                 y=1.02, fontsize=11.5, fontweight="bold")
    fig.tight_layout()
    save(fig, "Figure_3_effect_concordance")


def figure4_robustness() -> None:
    counts = pd.read_csv(OUTPUTS / "synthesis/differential_association_counts.csv")
    summary = pd.read_csv(OUTPUTS / "prediction/cv_summary.csv")
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.48, wspace=0.42)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    lcpm_components = ["qmp", "row_closed", "clr_minimum_positive", "clr_multiplicative"]
    lcpm_labels = ["QMP", "Row-closed", "CLR\nmin-pos.", "CLR\nmult."]
    pooled = counts.loc[
        (counts["cohort"] == "LCPM")
        & (counts["filter_specification"] == "pooled_outcome_blind")
    ].set_index("component")
    source = counts.loc[
        (counts["cohort"] == "LCPM")
        & (counts["filter_specification"] == "source_aligned_group_union")
    ].set_index("component")
    x = np.arange(len(lcpm_components))
    width = 0.36
    pvals = pooled.loc[lcpm_components, "significant_q_lt_0_05"].to_numpy()
    svals = source.loc[lcpm_components, "significant_q_lt_0_05"].to_numpy()
    ax1.bar(x - width / 2, pvals, width, color=COLORS["QMP"], label="Pooled outcome-blind")
    ax1.bar(x + width / 2, svals, width, color="#9CA3AF", label="Source-aligned group union")
    for offset, values in [(-width / 2, pvals), (width / 2, svals)]:
        for pos, value in zip(x + offset, values):
            ax1.text(pos, value + 0.25, str(int(value)), ha="center", fontsize=7.8)
    ax1.set_xticks(x, lcpm_labels, fontsize=7.0)
    ax1.set_ylim(0, 11.0)
    ax1.set_ylabel("Features with BH q < 0.05")
    ax1.set_title("LCPM filtering sensitivity")
    ax1.legend(frameon=False, fontsize=6.8, loc="upper center")
    ax1.grid(axis="y", color="#E5E7EB", linewidth=0.6)
    ax1.text(-0.14, 1.01, "A", transform=ax1.transAxes, fontsize=12, fontweight="bold")

    meta_components = [
        "prevalence", "qmp_nonzero", "row_closed_nonzero",
        "clr_minimum_positive", "clr_multiplicative",
    ]
    meta_labels = ["Prev.", "QMP", "Row-\nclosed", "CLR\nmin", "CLR\nmult"]
    core = counts.loc[
        (counts["cohort"] == "MetaCardis") & (counts["adjustment"] == "core")
    ].set_index("component")
    meds = counts.loc[
        (counts["cohort"] == "MetaCardis")
        & (counts["adjustment"] == "core_plus_medications")
    ].set_index("component")
    x2 = np.arange(len(meta_components))
    core_values = core.loc[meta_components, "significant_q_lt_0_05"].to_numpy()
    med_values = meds.loc[meta_components, "significant_q_lt_0_05"].to_numpy()
    ax2.bar(x2 - width / 2, core_values, width, color=COLORS["Clinical"], label="Core")
    ax2.bar(x2 + width / 2, med_values, width, color="#B8A6C9", label="Core + medications")
    for offset, values in [(-width / 2, core_values), (width / 2, med_values)]:
        for pos, value in zip(x2 + offset, values):
            ax2.text(pos, value + 0.35, str(int(value)), ha="center", fontsize=7.8)
    ax2.set_xticks(x2, meta_labels, fontsize=7.0)
    ax2.set_ylim(0, 18.0)
    ax2.set_title("MetaCardis covariate and CLR sensitivity")
    ax2.legend(frameon=False, fontsize=7.0, loc="upper center")
    ax2.grid(axis="y", color="#E5E7EB", linewidth=0.6)
    ax2.text(-0.20, 1.01, "B", transform=ax2.transAxes, fontsize=12, fontweight="bold")

    clinical_models = [
        "Clinical (HGB)", "QMP + clinical", "Row-closed + clinical", "CLR + clinical"
    ]
    display = ["Clinical only\n(HGB)", "QMP +\nclinical", "Row-closed +\nclinical", "CLR +\nclinical"]
    meta_summary = summary.loc[
        (summary["cohort"] == "MetaCardis") & summary["model"].isin(clinical_models)
    ].set_index("model")
    values = meta_summary.loc[clinical_models, "roc_auc_mean"].to_numpy()
    sds = meta_summary.loc[clinical_models, "roc_auc_sd"].to_numpy()
    bar_colors = [COLORS["Clinical"], COLORS["QMP"], COLORS["Row-closed"], COLORS["CLR"]]
    ax3.bar(range(4), values, yerr=sds, color=bar_colors, alpha=0.9, capsize=3)
    ax3.set_xticks(range(4), display, fontsize=8)
    ax3.set_ylim(0.86, 0.91)
    ax3.set_ylabel("Repeated-CV ROC AUC")
    ax3.set_title("Estimator-matched MetaCardis clinical comparison")
    ax3.grid(axis="y", color="#E5E7EB", linewidth=0.6)
    for pos, (value, sd) in enumerate(zip(values, sds)):
        ax3.text(pos, value + sd + 0.0012, f"{value:.3f}", ha="center", fontsize=8.2)
    ax3.text(-0.12, 1.01, "C", transform=ax3.transAxes, fontsize=12, fontweight="bold")

    fig.suptitle("Filtering, zero handling, and clinical structure materially qualify conclusions",
                 y=0.995, fontsize=11.3, fontweight="bold")
    save(fig, "Figure_4_robustness")


def supplementary_figures() -> None:
    shared = pd.read_csv(OUTPUTS / "synthesis/shared_exact_species.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    comparisons = [
        (
            "qmp_minus_row_closed_z_gap_lcpm",
            "qmp_minus_row_closed_z_gap_metacardis",
            "QMP − row-closed",
        ),
        ("qmp_minus_clr_z_gap_lcpm", "qmp_minus_clr_z_gap_metacardis", "QMP − CLR"),
    ]
    for ax, (xcol, ycol, title) in zip(axes, comparisons):
        r = float(shared[xcol].corr(shared[ycol]))
        ax.scatter(shared[xcol], shared[ycol], s=20, color=COLORS["QMP"], alpha=0.65,
                   edgecolor="white", linewidth=0.4)
        ax.axhline(0, color="#D1D5DB", linewidth=0.8)
        ax.axvline(0, color="#D1D5DB", linewidth=0.8)
        ax.set_xlabel("LCPM standardized representation gap")
        ax.set_ylabel("MetaCardis standardized representation gap")
        ax.set_title(f"{title}; Pearson r = {r:.3f}")
    panel_letter(axes[0], "A")
    panel_letter(axes[1], "B")
    fig.suptitle(f"Representation sensitivity was weakly correlated across {len(shared)} exact shared species",
                 y=1.02, fontsize=11.2, fontweight="bold")
    fig.tight_layout()
    save(fig, "Supplementary_Figure_1_shared_species")

    global_da = pd.read_csv(OUTPUTS / "associations/lcpm_global_prevalence_sensitivity.csv")
    summarized = (global_da.groupby(["prevalence_threshold", "representation", "test"], as_index=False)
                  .agg(tested_features=("tested_features", "first"), calls=("significant_q_lt_0_05", "sum")))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    for rep in ["QMP", "Row-closed"]:
        sub = summarized.loc[(summarized["representation"] == rep) & (summarized["test"] == "Kruskal-Wallis")]
        axes[0].plot(sub["prevalence_threshold"] * 100, sub["calls"], marker="o",
                     color=COLORS[rep], label=rep)
    for test, ls in [("one-way ANOVA", "-"), ("Kruskal-Wallis", "--")]:
        sub = summarized.loc[(summarized["representation"] == "CLR") & (summarized["test"] == test)]
        axes[0].plot(sub["prevalence_threshold"] * 100, sub["calls"], marker="o", linestyle=ls,
                     color=COLORS["CLR"], label=f"CLR ({'ANOVA' if test == 'one-way ANOVA' else 'KW'})")
    axes[0].set_xlabel("Prevalence threshold (%)")
    axes[0].set_ylabel("Features with BH q < 0.05")
    axes[0].set_title("LCPM three-group discovery sensitivity")
    axes[0].legend(frameon=False, fontsize=7.5)
    panel_letter(axes[0], "A")
    feature_counts = summarized.drop_duplicates("prevalence_threshold").sort_values("prevalence_threshold")
    axes[1].plot(feature_counts["prevalence_threshold"] * 100, feature_counts["tested_features"],
                 marker="o", color=COLORS["Muted"])
    axes[1].set_xlabel("Prevalence threshold (%)")
    axes[1].set_ylabel("Features tested")
    axes[1].set_title("Analysis-set size")
    panel_letter(axes[1], "B")
    fig.tight_layout()
    save(fig, "Supplementary_Figure_2_prevalence_sensitivity")

    load = pd.read_csv(ROOT / "paper/analysis/load_sensitivity_summary.csv")
    fig, ax = plt.subplots(figsize=(5.4, 2.7))
    y = np.arange(len(load))
    auc = load["direct_load_auc"].to_numpy()
    low = load["bootstrap_auc_ci95_low"].to_numpy()
    high = load["bootstrap_auc_ci95_high"].to_numpy()
    ax.errorbar(
        auc,
        y,
        xerr=np.vstack([auc - low, high - auc]),
        fmt="o",
        ms=6,
        color=COLORS["QMP"],
        ecolor=COLORS["QMP"],
        capsize=3,
        linewidth=1.4,
    )
    ax.axvline(0.5, color="#9CA3AF", linestyle="--", linewidth=0.9)
    ax.set_yticks(y, ["LCPM", "MetaCardis"])
    ax.set_xlim(0.35, 0.65)
    ax.set_xlabel("ROC AUC of total microbial load alone")
    ax.set_title("Post hoc microbial-load discrimination")
    for yy, value in zip(y, auc):
        ax.text(value, yy + 0.20, f"{value:.3f}", ha="center", fontsize=8)
    fig.tight_layout()
    save(fig, "Supplementary_Figure_3_total_load")


def main() -> None:
    style()
    figure1_workflow()
    figure2_prediction()
    figure3_effect_concordance()
    figure4_robustness()
    supplementary_figures()
    print(f"Created figures in {FIGURES}")


if __name__ == "__main__":
    main()
