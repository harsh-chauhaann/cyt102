#!/usr/bin/env python3
"""
render_binary_rankings.py

Top-vs-bottom contender analysis for binary M-L delta_E results.

Reads:
    mp_pipeline_results/delta_study_builders/all_deltaE_combined.csv

Writes:
    mp_pipeline_results/comparison_plots_top_bottom_binary/

The script focuses on both ends of the accepted delta_E distribution:
- Top contenders: most negative delta_E
- Bottom contenders: weakest accepted dips, i.e. delta_E values closest to 0
"""

from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# PATHS
# ============================================================
ROOT = Path(__file__).resolve().parent
INPUT_CSV = ROOT / "delta_analysis" / "all_deltaE_combined.csv"
OUT_DIR = ROOT / "comparison_plots_top_bottom_binary"

SUMMARY_DIR = OUT_DIR / "summary"
RANKING_DIR = OUT_DIR / "rankings"
SCATTER_DIR = OUT_DIR / "scatter"
CSV_DIR = OUT_DIR / "tables"

for directory in [OUT_DIR, SUMMARY_DIR, RANKING_DIR, SCATTER_DIR, CSV_DIR]:
    directory.mkdir(parents=True, exist_ok=True)


# ============================================================
# PARAMETERS
# ============================================================
TOP_N_CONFIGS = 20
TOP_N_CHEMSYS = 15
ANNOTATE_POINTS = 8


# ============================================================
# HELPERS
# ============================================================
def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))


def save_fig(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def apply_negative_xlim(ax, values, right_pad=0.02, pad_frac=0.08):
    vals = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if vals.empty:
        ax.set_xlim(-1.0, 0.02)
        return

    left = float(vals.min())
    right = max(float(vals.max()), 0.0)
    span = max(right - left, 1e-6)
    ax.set_xlim(left - pad_frac * span, right + right_pad)


def parse_binary_chemsys(chemsys: str):
    parts = str(chemsys).split("-")
    if len(parts) != 2:
        raise ValueError(f"Expected binary chemsys like 'Sc-O', got: {chemsys}")
    return parts[0], parts[1]


def add_descriptor_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    parsed = out["chemsys"].apply(parse_binary_chemsys)
    out["metal"] = parsed.str[0]
    out["ligand"] = parsed.str[1]
    return out


def top_configurations(df: pd.DataFrame, n: int = TOP_N_CONFIGS) -> pd.DataFrame:
    return df.sort_values("delta_E", ascending=True).head(n).copy()


def bottom_configurations(df: pd.DataFrame, n: int = TOP_N_CONFIGS) -> pd.DataFrame:
    return df.sort_values("delta_E", ascending=False).head(n).copy()


def chemsys_summary(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["chemsys", "metal", "ligand"], as_index=False)
          .agg(
              n=("delta_E", "size"),
              best_delta_E=("delta_E", "min"),
              worst_delta_E=("delta_E", "max"),
              mean_delta_E=("delta_E", "mean"),
              median_delta_E=("delta_E", "median"),
              std_delta_E=("delta_E", "std"),
          )
    )


# ============================================================
# TABLES
# ============================================================
def write_tables(df: pd.DataFrame):
    summary = chemsys_summary(df)
    summary.to_csv(CSV_DIR / "system_level_overview.csv", index=False)

    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)

    top_cfg.to_csv(CSV_DIR / "top_configurations.csv", index=False)
    bottom_cfg.to_csv(CSV_DIR / "bottom_configurations.csv", index=False)

    strongest = summary.sort_values("best_delta_E", ascending=True).head(TOP_N_CHEMSYS)
    weakest = summary.sort_values("best_delta_E", ascending=False).head(TOP_N_CHEMSYS)

    strongest.to_csv(CSV_DIR / "strongest_chemsys.csv", index=False)
    weakest.to_csv(CSV_DIR / "weakest_chemsys.csv", index=False)

    ligand_extremes = pd.DataFrame([
        {
            "group": "top_configurations",
            "ligand": ligand,
            "count": count,
        }
        for ligand, count in top_cfg["ligand"].value_counts().sort_index().items()
    ] + [
        {
            "group": "bottom_configurations",
            "ligand": ligand,
            "count": count,
        }
        for ligand, count in bottom_cfg["ligand"].value_counts().sort_index().items()
    ])
    ligand_extremes.to_csv(CSV_DIR / "ligand_extreme_counts.csv", index=False)

    metal_extremes = pd.DataFrame([
        {
            "group": "top_configurations",
            "metal": metal,
            "count": count,
        }
        for metal, count in top_cfg["metal"].value_counts().sort_index().items()
    ] + [
        {
            "group": "bottom_configurations",
            "metal": metal,
            "count": count,
        }
        for metal, count in bottom_cfg["metal"].value_counts().sort_index().items()
    ])
    metal_extremes.to_csv(CSV_DIR / "metal_extreme_counts.csv", index=False)

    with open(SUMMARY_DIR / "quick_summary.txt", "w", encoding="utf-8") as handle:
        handle.write(f"rows: {len(df)}\n")
        handle.write(f"unique chemsys: {df['chemsys'].nunique()}\n")
        handle.write(f"unique metals: {df['metal'].nunique()}\n")
        handle.write(f"ligands: {', '.join(sorted(df['ligand'].dropna().unique()))}\n")
        handle.write(f"delta_E min: {df['delta_E'].min():.6f}\n")
        handle.write(f"delta_E max: {df['delta_E'].max():.6f}\n")
        handle.write(f"delta_E mean: {df['delta_E'].mean():.6f}\n\n")

        best_row = top_cfg.iloc[0]
        worst_row = bottom_cfg.iloc[0]
        handle.write("Best overall configuration:\n")
        handle.write(best_row.to_string())
        handle.write("\n\n")
        handle.write("Weakest accepted configuration:\n")
        handle.write(worst_row.to_string())
        handle.write("\n\n")

        handle.write(f"Top {len(top_cfg)} configurations by most negative delta_E:\n")
        handle.write(
            top_cfg[["chemsys", "metal", "ligand", "source", "delta_E", "dip_cell_length", "baseline_slope"]]
            .to_string(index=False)
        )
        handle.write("\n\n")

        handle.write(f"Bottom {len(bottom_cfg)} configurations by weakest accepted delta_E:\n")
        handle.write(
            bottom_cfg[["chemsys", "metal", "ligand", "source", "delta_E", "dip_cell_length", "baseline_slope"]]
            .to_string(index=False)
        )
        handle.write("\n")


# ============================================================
# PLOTS
# ============================================================
def plot_configuration_extremes(df: pd.DataFrame):
    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)

    top_labels = [f"{row.chemsys} | {row.source}" for _, row in top_cfg.iterrows()]
    bottom_labels = [f"{row.chemsys} | {row.source}" for _, row in bottom_cfg.iterrows()]

    fig_h = max(7, 0.30 * max(len(top_cfg), len(bottom_cfg)))
    fig, axes = plt.subplots(1, 2, figsize=(18, fig_h), sharex=False)

    top_vals = top_cfg["delta_E"].values
    top_y = np.arange(len(top_cfg))
    axes[0].barh(top_y, top_vals, color="#1f77b4")
    axes[0].set_yticks(top_y)
    axes[0].set_yticklabels(top_labels, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("delta_E")
    axes[0].set_title(f"Top {len(top_cfg)} contenders")
    apply_negative_xlim(axes[0], top_vals)

    bottom_vals = bottom_cfg["delta_E"].values
    bottom_y = np.arange(len(bottom_cfg))
    axes[1].barh(bottom_y, bottom_vals, color="#d62728")
    axes[1].set_yticks(bottom_y)
    axes[1].set_yticklabels(bottom_labels, fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("delta_E")
    axes[1].set_title(f"Bottom {len(bottom_cfg)} accepted contenders")
    apply_negative_xlim(axes[1], bottom_vals)

    save_fig(fig, RANKING_DIR / "configuration_extremes_top_vs_bottom.png")


def plot_chemsys_extremes(df: pd.DataFrame):
    summary = chemsys_summary(df)
    strongest = summary.sort_values("best_delta_E", ascending=True).head(TOP_N_CHEMSYS)
    weakest = summary.sort_values("best_delta_E", ascending=False).head(TOP_N_CHEMSYS)

    fig_h = max(6, 0.34 * max(len(strongest), len(weakest)))
    fig, axes = plt.subplots(1, 2, figsize=(16, fig_h), sharex=False)

    left_vals = strongest["best_delta_E"].values
    left_y = np.arange(len(strongest))
    axes[0].barh(left_y, left_vals, color="#2ca02c")
    axes[0].set_yticks(left_y)
    axes[0].set_yticklabels(strongest["chemsys"], fontsize=9)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Best delta_E")
    axes[0].set_title(f"Strongest {len(strongest)} M-L systems")
    apply_negative_xlim(axes[0], left_vals)

    right_vals = weakest["best_delta_E"].values
    right_y = np.arange(len(weakest))
    axes[1].barh(right_y, right_vals, color="#ff7f0e")
    axes[1].set_yticks(right_y)
    axes[1].set_yticklabels(weakest["chemsys"], fontsize=9)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Best delta_E")
    axes[1].set_title(f"Weakest {len(weakest)} M-L systems")
    apply_negative_xlim(axes[1], right_vals)

    save_fig(fig, RANKING_DIR / "chemsys_extremes_top_vs_bottom.png")


def plot_overview_scatter(df: pd.DataFrame):
    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)

    fig, ax = plt.subplots(figsize=(11, 7))
    ax.scatter(
        df["dip_cell_length"],
        df["delta_E"],
        s=28,
        color="lightgray",
        alpha=0.65,
        label="All accepted configurations",
    )
    ax.scatter(
        top_cfg["dip_cell_length"],
        top_cfg["delta_E"],
        s=70,
        color="#1f77b4",
        label="Top contenders",
    )
    ax.scatter(
        bottom_cfg["dip_cell_length"],
        bottom_cfg["delta_E"],
        s=70,
        color="#d62728",
        label="Bottom contenders",
    )

    for _, row in top_cfg.head(min(ANNOTATE_POINTS, len(top_cfg))).iterrows():
        ax.annotate(
            row["chemsys"],
            (row["dip_cell_length"], row["delta_E"]),
            fontsize=8,
            xytext=(4, 3),
            textcoords="offset points",
        )

    for _, row in bottom_cfg.head(min(ANNOTATE_POINTS, len(bottom_cfg))).iterrows():
        ax.annotate(
            row["chemsys"],
            (row["dip_cell_length"], row["delta_E"]),
            fontsize=8,
            xytext=(4, -10),
            textcoords="offset points",
        )

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Dip cell length")
    ax.set_ylabel("delta_E")
    ax.set_title("Accepted M-L configurations: strongest and weakest contenders")
    ax.legend()

    y_vals = pd.concat([df["delta_E"], pd.Series([0.0])], ignore_index=True)
    vals = pd.to_numeric(y_vals, errors="coerce").dropna()
    ymin = float(vals.min()) if not vals.empty else -1.0
    ymax = float(vals.max()) if not vals.empty else 0.02
    pad = max(0.03, 0.08 * max(abs(ymin), abs(ymax), 1e-6))
    ax.set_ylim(ymin - pad, ymax + pad)

    save_fig(fig, SCATTER_DIR / "delta_vs_dip_cell_length_top_bottom.png")


def plot_ligand_extreme_counts(df: pd.DataFrame):
    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)

    ligands = sorted(df["ligand"].dropna().unique())
    top_counts = [int((top_cfg["ligand"] == ligand).sum()) for ligand in ligands]
    bottom_counts = [int((bottom_cfg["ligand"] == ligand).sum()) for ligand in ligands]

    x = np.arange(len(ligands))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - width / 2, top_counts, width=width, color="#1f77b4", label="Top")
    ax.bar(x + width / 2, bottom_counts, width=width, color="#d62728", label="Bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(ligands)
    ax.set_ylabel("Count")
    ax.set_title("Ligand makeup of top vs bottom contenders")
    ax.legend()

    save_fig(fig, SUMMARY_DIR / "ligand_top_bottom_counts.png")


def plot_metal_extreme_counts(df: pd.DataFrame):
    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)

    top_counts = top_cfg["metal"].value_counts().head(10).sort_values(ascending=True)
    bottom_counts = bottom_cfg["metal"].value_counts().head(10).sort_values(ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=False)

    axes[0].barh(top_counts.index, top_counts.values, color="#1f77b4")
    axes[0].set_xlabel("Count")
    axes[0].set_title("Most frequent metals in top contenders")

    axes[1].barh(bottom_counts.index, bottom_counts.values, color="#d62728")
    axes[1].set_xlabel("Count")
    axes[1].set_title("Most frequent metals in bottom contenders")

    save_fig(fig, SUMMARY_DIR / "metal_top_bottom_counts.png")


# ============================================================
# MAIN
# ============================================================
def main():
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input CSV not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)
    if df.empty:
        raise ValueError("Input CSV is empty.")

    required = {"chemsys", "source", "delta_E", "dip_cell_length", "baseline_slope"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = add_descriptor_columns(df)

    numeric_cols = ["delta_E", "dip_cell_length", "baseline_slope"]
    for column in numeric_cols:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(subset=["delta_E", "dip_cell_length"]).copy()
    if df.empty:
        raise ValueError("No rows remain after numeric cleanup.")

    print(f"Loaded rows: {len(df)}")
    print(f"Unique chemsys: {df['chemsys'].nunique()}")
    print(f"Unique metals: {df['metal'].nunique()}")
    print(f"Ligands: {sorted(df['ligand'].dropna().unique())}")

    write_tables(df)
    plot_configuration_extremes(df)
    plot_chemsys_extremes(df)
    plot_overview_scatter(df)
    plot_ligand_extreme_counts(df)
    plot_metal_extreme_counts(df)

    print("\nDone.")
    print(f"Saved top-vs-bottom analysis to: {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()

