#!/usr/bin/env python3
"""Generate publication-quality figures for the cancer driver antibody design preprint.

Reads from: results/analysis/all_scores.csv, results/analysis/cross_campaign_summary.json
Outputs to:  out/figures/

Uses Paul Tol color-blind friendly palette per Saloni Dattani viz guidelines.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

OUT = Path("out/figures")
OUT.mkdir(parents=True, exist_ok=True)

# Paul Tol vibrant palette (colorblind-safe)
TOLVIB = {
    "blue":    "#0077BB",
    "cyan":    "#33BBEE",
    "teal":    "#009988",
    "orange":  "#EE7733",
    "red":     "#CC3311",
    "magenta": "#EE3377",
    "grey":    "#BBBBBB",
}
PAL10 = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
    "#AA3377", "#BBBBBB", "#EE7733", "#0077BB", "#009988",
]

CAMPAIGN_LABELS = {
    "b7h3_vhh": "B7-H3",
    "cd47_vhh": "CD47",
    "ceacam5_vhh": "CEACAM5",
    "egfr_cetuximab_vhh": "EGFR",
    "egfrviii_vhh": "EGFRvIII",
    "epha2_vhh": "EphA2",
    "gpc2_vhh": "GPC2",
    "her2_domIV_vhh": "HER2-DIV",
    "msln_cterm_vhh": "MSLN-C",
    "msln_nterm_vhh": "MSLN-N",
}
CAMPAIGN_ORDER = list(CAMPAIGN_LABELS.keys())

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 200,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load_data():
    df = pd.read_csv("results/analysis/all_scores.csv")
    with open("results/analysis/cross_campaign_summary.json") as f:
        summary = json.load(f)
    return df, {s["campaign"]: s for s in summary}


def fig1_pass_rate_bar(summary):
    """Horizontal bar chart of filter pass rates by campaign."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    campaigns = CAMPAIGN_ORDER
    labels = [CAMPAIGN_LABELS[c] for c in campaigns]
    rates = [summary[c]["pass_rate_pct"] for c in campaigns]
    counts = [summary[c]["passed_filters"] for c in campaigns]

    colors = [PAL10[i] for i in range(len(campaigns))]
    bars = ax.barh(range(len(campaigns)), rates, color=colors, edgecolor="white", linewidth=0.5)

    for i, (bar, count, rate) in enumerate(zip(bars, counts, rates)):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f"{count} ({rate:.1f}%)", va="center", fontsize=8.5)

    ax.set_yticks(range(len(campaigns)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Filter pass rate (%)")
    ax.set_title("Design Filter Pass Rates Across Cancer Targets")
    ax.invert_yaxis()
    ax.set_xlim(0, max(rates) * 1.35)
    fig.savefig(OUT / "fig1_pass_rates.pdf")
    fig.savefig(OUT / "fig1_pass_rates.png")
    plt.close(fig)
    print("  fig1_pass_rates")


def fig2_pae_distributions(df):
    """Box + strip plot of pae distributions per campaign."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    data = []
    labels = []
    for i, c in enumerate(CAMPAIGN_ORDER):
        vals = df[df["campaign"] == c]["pae"].dropna()
        data.append(vals.values)
        labels.append(CAMPAIGN_LABELS[c])

    bp = ax.boxplot(data, vert=True, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", linewidth=1.5),
                    flierprops=dict(marker=".", markersize=2, alpha=0.3))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(PAL10[i])
        patch.set_alpha(0.7)

    ax.axhline(y=10, color=TOLVIB["red"], linestyle="--", linewidth=1, alpha=0.7, label="Filter threshold (pAE < 10)")
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("pAE (predicted aligned error)")
    ax.set_title("pAE Score Distributions by Target")
    ax.legend(loc="upper left", framealpha=0.9)
    fig.savefig(OUT / "fig2_pae_distributions.pdf")
    fig.savefig(OUT / "fig2_pae_distributions.png")
    plt.close(fig)
    print("  fig2_pae_distributions")


def fig3_rmsd_distributions(df):
    """Box + strip plot of target-aligned CDR RMSD per campaign."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    data = []
    labels = []
    for c in CAMPAIGN_ORDER:
        vals = df[df["campaign"] == c]["rmsd"].dropna()
        data.append(vals.values)
        labels.append(CAMPAIGN_LABELS[c])

    bp = ax.boxplot(data, vert=True, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", linewidth=1.5),
                    flierprops=dict(marker=".", markersize=2, alpha=0.3))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(PAL10[i])
        patch.set_alpha(0.7)

    ax.axhline(y=2.0, color=TOLVIB["red"], linestyle="--", linewidth=1, alpha=0.7, label="Filter threshold (RMSD < 2.0 Å)")
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Target-aligned CDR RMSD (Å)")
    ax.set_title("CDR RMSD Distributions by Target")
    ax.legend(loc="upper left", framealpha=0.9)
    fig.savefig(OUT / "fig3_rmsd_distributions.pdf")
    fig.savefig(OUT / "fig3_rmsd_distributions.png")
    plt.close(fig)
    print("  fig3_rmsd_distributions")


def fig4_pae_vs_rmsd_scatter(df):
    """Joint scatter of pae vs rmsd colored by campaign with quadrant annotations."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, c in enumerate(CAMPAIGN_ORDER):
        sub = df[df["campaign"] == c]
        ax.scatter(sub["pae"], sub["rmsd"], c=PAL10[i], label=CAMPAIGN_LABELS[c],
                   s=8, alpha=0.4, edgecolors="none")

    ax.axvline(x=10, color=TOLVIB["red"], linestyle="--", linewidth=1, alpha=0.6)
    ax.axhline(y=2.0, color=TOLVIB["red"], linestyle="--", linewidth=1, alpha=0.6)

    ax.fill_between([0, 10], 0, 2.0, alpha=0.06, color=TOLVIB["teal"])
    ax.text(5, 1.0, "PASS", fontsize=14, fontweight="bold", color=TOLVIB["teal"],
            ha="center", va="center", alpha=0.5)

    ax.set_xlabel("pAE (predicted aligned error)")
    ax.set_ylabel("Target-aligned CDR RMSD (Å)")
    ax.set_title("Joint pAE–RMSD Distribution (all 4,085 designs)")
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.legend(loc="upper right", markerscale=3, framealpha=0.9, fontsize=8, ncol=2)
    fig.savefig(OUT / "fig4_pae_vs_rmsd.pdf")
    fig.savefig(OUT / "fig4_pae_vs_rmsd.png")
    plt.close(fig)
    print("  fig4_pae_vs_rmsd")


def fig5_design_funnel(summary):
    """Funnel diagram showing attrition from backbones → sequences → passed."""
    fig, ax = plt.subplots(figsize=(7, 5))

    total_backbones = 817
    total_sequences = 4085
    total_passed = sum(s["passed_filters"] for s in summary.values())

    stages = ["Backbones\n(Stage 1: RFdiffusion)", "Sequences\n(Stage 2: ProteinMPNN)", "Scored & Filtered\n(Stage 3: RF2)"]
    values = [total_backbones, total_sequences, total_passed]
    colors_f = [TOLVIB["blue"], TOLVIB["cyan"], TOLVIB["teal"]]

    max_w = 0.8
    widths = [max_w * v / max(values) for v in values]

    for i, (stage, val, w, col) in enumerate(zip(stages, values, widths, colors_f)):
        y = 2 - i
        rect = FancyBboxPatch(((1-w)/2, y-0.35), w, 0.7, boxstyle="round,pad=0.05",
                              facecolor=col, edgecolor="white", linewidth=2, alpha=0.85)
        ax.add_patch(rect)
        ax.text(0.5, y, f"{val:,}", ha="center", va="center",
                fontsize=18, fontweight="bold", color="white")
        ax.text(0.5, y - 0.5, stage, ha="center", va="top", fontsize=9, color="#333")

        if i < len(values) - 1:
            ratio = values[i+1] / values[i]
            ax.annotate(f"×{ratio:.1f}", xy=(0.92, y - 0.35), fontsize=8, color="#666")

    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-1, 3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Design Funnel: 10 Cancer Targets", fontsize=13, fontweight="bold", pad=15)
    fig.savefig(OUT / "fig5_design_funnel.pdf")
    fig.savefig(OUT / "fig5_design_funnel.png")
    plt.close(fig)
    print("  fig5_design_funnel")


def fig6_metric_heatmap(summary):
    """Heatmap of key metrics across campaigns."""
    campaigns = CAMPAIGN_ORDER
    labels = [CAMPAIGN_LABELS[c] for c in campaigns]
    metrics = ["pae_median", "rmsd_median", "interaction_pae_median", "pred_lddt_median",
               "framework_aligned_H3_rmsd_median"]
    metric_labels = ["pAE\n(median)", "CDR RMSD\n(median, Å)", "iPAE\n(median)",
                     "pLDDT\n(median)", "H3 RMSD\n(median, Å)"]

    data = np.zeros((len(campaigns), len(metrics)))
    for i, c in enumerate(campaigns):
        for j, m in enumerate(metrics):
            data[i, j] = summary[c].get(m, np.nan)

    fig, ax = plt.subplots(figsize=(8, 5))

    norm_data = np.zeros_like(data)
    for j in range(data.shape[1]):
        col = data[:, j]
        lo, hi = np.nanmin(col), np.nanmax(col)
        if hi - lo > 1e-12:
            norm_data[:, j] = (col - lo) / (hi - lo)

    # For pLDDT, invert (higher is better)
    if 3 < data.shape[1]:
        norm_data[:, 3] = 1 - norm_data[:, 3]

    im = ax.imshow(norm_data, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=1)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if norm_data[i,j] > 0.65 else "black")

    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metric_labels, fontsize=9)
    ax.set_yticks(range(len(campaigns)))
    ax.set_yticklabels(labels)
    ax.set_title("Structural Quality Metrics Across Targets\n(green = better, red = worse)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Relative score (lower = better)")
    fig.savefig(OUT / "fig6_metric_heatmap.pdf")
    fig.savefig(OUT / "fig6_metric_heatmap.png")
    plt.close(fig)
    print("  fig6_metric_heatmap")


def fig7_h3_rmsd_violin(df):
    """Violin plot of framework-aligned H3 RMSD — key for CDR-H3 loop accuracy."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    data = []
    positions = []
    for i, c in enumerate(CAMPAIGN_ORDER):
        vals = df[df["campaign"] == c]["framework_aligned_H3_rmsd"].dropna()
        if len(vals) > 1:
            data.append(vals.values)
            positions.append(i)

    if data:
        vp = ax.violinplot(data, positions=positions, showmedians=True, widths=0.7)
        for i, body in enumerate(vp["bodies"]):
            body.set_facecolor(PAL10[positions[i]])
            body.set_alpha(0.7)
        vp["cmedians"].set_color("black")

    ax.set_xticks(range(len(CAMPAIGN_ORDER)))
    ax.set_xticklabels([CAMPAIGN_LABELS[c] for c in CAMPAIGN_ORDER], rotation=45, ha="right")
    ax.set_ylabel("Framework-aligned H3 RMSD (Å)")
    ax.set_title("CDR-H3 Loop Accuracy by Target")
    fig.savefig(OUT / "fig7_h3_rmsd_violin.pdf")
    fig.savefig(OUT / "fig7_h3_rmsd_violin.png")
    plt.close(fig)
    print("  fig7_h3_rmsd_violin")


def fig8_lddt_distributions(df):
    """Box plot of predicted LDDT (global fold quality)."""
    fig, ax = plt.subplots(figsize=(8, 4))
    data = []
    for c in CAMPAIGN_ORDER:
        vals = df[df["campaign"] == c]["pred_lddt"].dropna()
        data.append(vals.values)

    bp = ax.boxplot(data, vert=True, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", linewidth=1.5),
                    flierprops=dict(marker=".", markersize=2, alpha=0.3))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(PAL10[i])
        patch.set_alpha(0.7)

    ax.set_xticklabels([CAMPAIGN_LABELS[c] for c in CAMPAIGN_ORDER], rotation=45, ha="right")
    ax.set_ylabel("Predicted LDDT")
    ax.set_title("Global Fold Quality (pLDDT) by Target")
    fig.savefig(OUT / "fig8_lddt_distributions.pdf")
    fig.savefig(OUT / "fig8_lddt_distributions.png")
    plt.close(fig)
    print("  fig8_lddt_distributions")


def fig9_top_candidates_table(summary):
    """Visual table of top candidate per campaign."""
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis("off")

    cols = ["Target", "Best Design", "pAE", "CDR RMSD", "iPAE", "pLDDT", "H3 RMSD", "Score"]
    rows = []
    for c in CAMPAIGN_ORDER:
        s = summary[c]
        top = s.get("top10", [{}])
        if top:
            t = top[0]
            rows.append([
                CAMPAIGN_LABELS[c],
                t.get("tag", "")[:30],
                f"{t.get('pae', 0):.2f}",
                f"{t.get('rmsd', 0):.2f}",
                f"{t.get('interaction_pae', 0):.2f}",
                f"{t.get('pred_lddt', 0):.2f}",
                f"{t.get('framework_aligned_H3_rmsd', 0):.2f}",
                f"{t.get('composite_score', 0):.3f}",
            ])

    table = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)

    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor(TOLVIB["blue"])
            cell.set_text_props(color="white", fontweight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#f0f0f0")

    ax.set_title("Top Candidate per Target (ranked by composite score)", fontsize=12, fontweight="bold", pad=20)
    fig.savefig(OUT / "fig9_top_candidates.pdf")
    fig.savefig(OUT / "fig9_top_candidates.png")
    plt.close(fig)
    print("  fig9_top_candidates")


def fig10_pass_by_backbone(df, summary):
    """Scatter: designs per backbone vs pass rate — shows which backbone topologies succeed."""
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, c in enumerate(CAMPAIGN_ORDER):
        s = summary[c]
        n_designs = s["total_designs"]
        n_passed = s["passed_filters"]
        n_backbones = n_designs // 5  # 5 sequences per backbone
        ax.scatter(n_backbones, s["pass_rate_pct"], c=PAL10[i], s=80, zorder=5,
                   edgecolors="white", linewidth=0.5)
        ax.annotate(CAMPAIGN_LABELS[c], (n_backbones, s["pass_rate_pct"]),
                    textcoords="offset points", xytext=(8, 3), fontsize=8)

    ax.set_xlabel("Number of backbone designs")
    ax.set_ylabel("Filter pass rate (%)")
    ax.set_title("Backbone Count vs. Design Success Rate")
    fig.savefig(OUT / "fig10_backbone_vs_passrate.pdf")
    fig.savefig(OUT / "fig10_backbone_vs_passrate.png")
    plt.close(fig)
    print("  fig10_backbone_vs_passrate")


def main():
    print("Loading data...")
    df, summary = load_data()
    print(f"  {len(df)} designs, {len(summary)} campaigns")

    print("Generating figures...")
    fig1_pass_rate_bar(summary)
    fig2_pae_distributions(df)
    fig3_rmsd_distributions(df)
    fig4_pae_vs_rmsd_scatter(df)
    fig5_design_funnel(summary)
    fig6_metric_heatmap(summary)
    fig7_h3_rmsd_violin(df)
    fig8_lddt_distributions(df)
    fig9_top_candidates_table(summary)
    fig10_pass_by_backbone(df, summary)

    print(f"\nAll 10 figures saved to {OUT}/")


if __name__ == "__main__":
    main()
