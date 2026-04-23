#!/usr/bin/env python3
"""
figure_ternary.py

Top-vs-bottom artifact analysis for ternary M1-M2-L delta_E results.

Reads:
    mp_pipeline_results/delta_analysis_ternary/all_ternary_deltaE_combined.csv

Writes:
    mp_pipeline_results/comparison_plots_top_bottom_ternary/

Interpretation:
- More negative delta_E means a larger suspicious D3-induced artifact
- delta_E values close to 0 mean weaker accepted artifact signatures
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
INPUT_CSV = ROOT / "delta_analysis_ternary" / "all_ternary_deltaE_combined.csv"
OUT_DIR = ROOT / "comparison_plots_top_bottom_ternary"

SUMMARY_DIR = OUT_DIR / "summary"
RANKING_DIR = OUT_DIR / "rankings"
PAIR_DIR = OUT_DIR / "pair_comparison"
SCATTER_DIR = OUT_DIR / "scatter"
HEATMAP_DIR = OUT_DIR / "heatmaps"
CSV_DIR = OUT_DIR / "tables"

for directory in [OUT_DIR, SUMMARY_DIR, RANKING_DIR, PAIR_DIR, SCATTER_DIR, HEATMAP_DIR, CSV_DIR]:
    directory.mkdir(parents=True, exist_ok=True)


# ============================================================
# PARAMETERS
# ============================================================
TOP_N_CONFIGS = 20
TOP_N_CHEMSYS = 15
TOP_N_PAIRS = 15
ANNOTATE_POINTS = 8

METAL_ORDER = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"
]

IONIZATION_ENERGY = {
    "Sc": 633.1, "Ti": 658.8, "V": 650.9, "Cr": 652.9, "Mn": 717.3,
    "Fe": 762.5, "Co": 760.4, "Ni": 737.1, "Cu": 745.5, "Zn": 906.4,
    "Y": 600.0, "Zr": 640.1, "Nb": 652.1, "Mo": 684.3, "Tc": 702.0,
    "Ru": 710.2, "Rh": 719.7, "Pd": 804.4, "Ag": 731.0, "Cd": 867.8,
    "Hf": 658.5, "Ta": 761.0, "W": 770.0, "Re": 760.0, "Os": 840.0,
    "Ir": 880.0, "Pt": 870.0, "Au": 890.1, "Hg": 1007.1,
}


# ============================================================
# HELPERS
# ============================================================
def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))


def canonical_pair(m1: str, m2: str) -> str:
    a, b = sorted([str(m1), str(m2)], key=lambda x: x.upper())
    return f"{a}-{b}"


def save_fig(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def apply_delta_xlim(ax, values, right_pad=0.02, pad_frac=0.08):
    vals = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if vals.empty:
        ax.set_xlim(-1.0, 0.02)
        return

    left = float(vals.min())
    right = max(float(vals.max()), 0.0)
    span = max(right - left, 1e-6)
    ax.set_xlim(left - pad_frac * span, right + right_pad)


def add_descriptor_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    if "pair" not in out.columns:
        out["pair"] = out.apply(lambda r: canonical_pair(r["metal_1"], r["metal_2"]), axis=1)

    if "IE_1" not in out.columns:
        out["IE_1"] = out["metal_1"].map(IONIZATION_ENERGY)
    if "IE_2" not in out.columns:
        out["IE_2"] = out["metal_2"].map(IONIZATION_ENERGY)
    if "IE_avg" not in out.columns:
        out["IE_avg"] = (out["IE_1"] + out["IE_2"]) / 2.0
    if "IE_diff" not in out.columns:
        out["IE_diff"] = (out["IE_1"] - out["IE_2"]).abs()

    return out


def top_configurations(df: pd.DataFrame, n: int = TOP_N_CONFIGS) -> pd.DataFrame:
    return df.sort_values("delta_E", ascending=True).head(n).copy()


def bottom_configurations(df: pd.DataFrame, n: int = TOP_N_CONFIGS) -> pd.DataFrame:
    return df.sort_values("delta_E", ascending=False).head(n).copy()


def chemsys_summary(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["chemsys", "ligand"], as_index=False)
          .agg(
              n=("delta_E", "size"),
              best_delta_E=("delta_E", "min"),
              worst_delta_E=("delta_E", "max"),
              mean_delta_E=("delta_E", "mean"),
              median_delta_E=("delta_E", "median"),
              std_delta_E=("delta_E", "std"),
          )
    )


def pair_summary(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["pair", "ligand"], as_index=False)
          .agg(
              n=("delta_E", "size"),
              best_delta_E=("delta_E", "min"),
              worst_delta_E=("delta_E", "max"),
              mean_delta_E=("delta_E", "mean"),
              median_delta_E=("delta_E", "median"),
              std_delta_E=("delta_E", "std"),
          )
    )


def strongest_chemsys(summary: pd.DataFrame, n: int = TOP_N_CHEMSYS) -> pd.DataFrame:
    return (
        summary.groupby("chemsys", as_index=False)
               .agg(
                   ligand_examples=("ligand", lambda s: ",".join(sorted(set(map(str, s))))),
                   n_total=("n", "sum"),
                   best_overall_delta_E=("best_delta_E", "min"),
                   worst_overall_delta_E=("worst_delta_E", "max"),
                   mean_overall_delta_E=("mean_delta_E", "mean"),
               )
               .sort_values("best_overall_delta_E", ascending=True)
               .head(n)
               .copy()
    )


def weakest_chemsys(summary: pd.DataFrame, n: int = TOP_N_CHEMSYS) -> pd.DataFrame:
    return (
        summary.groupby("chemsys", as_index=False)
               .agg(
                   ligand_examples=("ligand", lambda s: ",".join(sorted(set(map(str, s))))),
                   n_total=("n", "sum"),
                   best_overall_delta_E=("best_delta_E", "min"),
                   worst_overall_delta_E=("worst_delta_E", "max"),
                   mean_overall_delta_E=("mean_delta_E", "mean"),
               )
               .sort_values("best_overall_delta_E", ascending=False)
               .head(n)
               .copy()
    )


def strongest_pairs(summary: pd.DataFrame, n: int = TOP_N_PAIRS) -> pd.DataFrame:
    return (
        summary.groupby("pair", as_index=False)
               .agg(
                   ligand_examples=("ligand", lambda s: ",".join(sorted(set(map(str, s))))),
                   n_total=("n", "sum"),
                   best_overall_delta_E=("best_delta_E", "min"),
                   worst_overall_delta_E=("worst_delta_E", "max"),
                   mean_overall_delta_E=("mean_delta_E", "mean"),
               )
               .sort_values("best_overall_delta_E", ascending=True)
               .head(n)
               .copy()
    )


def weakest_pairs(summary: pd.DataFrame, n: int = TOP_N_PAIRS) -> pd.DataFrame:
    return (
        summary.groupby("pair", as_index=False)
               .agg(
                   ligand_examples=("ligand", lambda s: ",".join(sorted(set(map(str, s))))),
                   n_total=("n", "sum"),
                   best_overall_delta_E=("best_delta_E", "min"),
                   worst_overall_delta_E=("worst_delta_E", "max"),
                   mean_overall_delta_E=("mean_delta_E", "mean"),
               )
               .sort_values("best_overall_delta_E", ascending=False)
               .head(n)
               .copy()
    )


# ============================================================
# TABLES
# ============================================================
def write_tables(df: pd.DataFrame):
    cs = chemsys_summary(df)
    ps = pair_summary(df)
    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)
    strong_cs = strongest_chemsys(cs, TOP_N_CHEMSYS)
    weak_cs = weakest_chemsys(cs, TOP_N_CHEMSYS)
    strong_pairs = strongest_pairs(ps, TOP_N_PAIRS)
    weak_pairs = weakest_pairs(ps, TOP_N_PAIRS)

    cs.to_csv(CSV_DIR / "system_level_overview.csv", index=False)
    ps.to_csv(CSV_DIR / "pair_summary.csv", index=False)
    top_cfg.to_csv(CSV_DIR / "top_configurations.csv", index=False)
    bottom_cfg.to_csv(CSV_DIR / "bottom_configurations.csv", index=False)
    strong_cs.to_csv(CSV_DIR / "strongest_chemsys.csv", index=False)
    weak_cs.to_csv(CSV_DIR / "weakest_chemsys.csv", index=False)
    strong_pairs.to_csv(CSV_DIR / "strongest_pairs.csv", index=False)
    weak_pairs.to_csv(CSV_DIR / "weakest_pairs.csv", index=False)

    ligand_counts = pd.DataFrame([
        {"group": "top_configurations", "ligand": ligand, "count": count}
        for ligand, count in top_cfg["ligand"].value_counts().sort_index().items()
    ] + [
        {"group": "bottom_configurations", "ligand": ligand, "count": count}
        for ligand, count in bottom_cfg["ligand"].value_counts().sort_index().items()
    ])
    ligand_counts.to_csv(CSV_DIR / "ligand_extreme_counts.csv", index=False)

    with open(SUMMARY_DIR / "quick_summary.txt", "w", encoding="utf-8") as handle:
        handle.write(f"rows: {len(df)}\n")
        handle.write(f"unique chemsys: {df['chemsys'].nunique()}\n")
        handle.write(f"unique pairs: {df['pair'].nunique()}\n")
        handle.write(f"ligands: {', '.join(sorted(df['ligand'].dropna().unique()))}\n")
        handle.write(f"delta_E min: {df['delta_E'].min():.6f}\n")
        handle.write(f"delta_E max: {df['delta_E'].max():.6f}\n")
        handle.write(f"delta_E mean: {df['delta_E'].mean():.6f}\n\n")

        handle.write("Largest artifact configuration:\n")
        handle.write(top_cfg.iloc[0].to_string())
        handle.write("\n\n")

        handle.write("Smallest accepted artifact configuration:\n")
        handle.write(bottom_cfg.iloc[0].to_string())
        handle.write("\n\n")

        handle.write(f"Top {len(top_cfg)} configurations by most negative delta_E:\n")
        handle.write(
            top_cfg[["chemsys", "pair", "ligand", "source", "delta_E", "dip_cell_length", "IE_avg", "IE_diff"]]
                  .to_string(index=False)
        )
        handle.write("\n\n")

        handle.write(f"Bottom {len(bottom_cfg)} configurations by weakest accepted delta_E:\n")
        handle.write(
            bottom_cfg[["chemsys", "pair", "ligand", "source", "delta_E", "dip_cell_length", "IE_avg", "IE_diff"]]
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
    fig, axes = plt.subplots(1, 2, figsize=(20, fig_h), sharex=False)

    top_vals = top_cfg["delta_E"].values
    top_y = np.arange(len(top_cfg))
    axes[0].barh(top_y, top_vals, color="#1f77b4")
    axes[0].set_yticks(top_y)
    axes[0].set_yticklabels(top_labels, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("delta_E")
    axes[0].set_title(f"Largest {len(top_cfg)} artifact configurations")
    apply_delta_xlim(axes[0], top_vals)

    bottom_vals = bottom_cfg["delta_E"].values
    bottom_y = np.arange(len(bottom_cfg))
    axes[1].barh(bottom_y, bottom_vals, color="#d62728")
    axes[1].set_yticks(bottom_y)
    axes[1].set_yticklabels(bottom_labels, fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("delta_E")
    axes[1].set_title(f"Smallest {len(bottom_cfg)} accepted artifact configurations")
    apply_delta_xlim(axes[1], bottom_vals)

    save_fig(fig, RANKING_DIR / "configuration_extremes_top_vs_bottom.png")


def plot_chemsys_extremes(df: pd.DataFrame):
    summary = chemsys_summary(df)
    strong = strongest_chemsys(summary, TOP_N_CHEMSYS)
    weak = weakest_chemsys(summary, TOP_N_CHEMSYS)

    fig_h = max(6, 0.34 * max(len(strong), len(weak)))
    fig, axes = plt.subplots(1, 2, figsize=(18, fig_h), sharex=False)

    left_vals = strong["best_overall_delta_E"].values
    left_y = np.arange(len(strong))
    axes[0].barh(left_y, left_vals, color="#2ca02c")
    axes[0].set_yticks(left_y)
    axes[0].set_yticklabels(strong["chemsys"], fontsize=9)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Best delta_E")
    axes[0].set_title(f"Largest-artifact {len(strong)} chemsys")
    apply_delta_xlim(axes[0], left_vals)

    right_vals = weak["best_overall_delta_E"].values
    right_y = np.arange(len(weak))
    axes[1].barh(right_y, right_vals, color="#ff7f0e")
    axes[1].set_yticks(right_y)
    axes[1].set_yticklabels(weak["chemsys"], fontsize=9)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Best delta_E")
    axes[1].set_title(f"Smallest-artifact {len(weak)} chemsys")
    apply_delta_xlim(axes[1], right_vals)

    save_fig(fig, RANKING_DIR / "chemsys_extremes_top_vs_bottom.png")


def plot_pair_extremes(df: pd.DataFrame):
    summary = pair_summary(df)
    strong = strongest_pairs(summary, TOP_N_PAIRS)
    weak = weakest_pairs(summary, TOP_N_PAIRS)

    fig_h = max(6, 0.34 * max(len(strong), len(weak)))
    fig, axes = plt.subplots(1, 2, figsize=(18, fig_h), sharex=False)

    left_vals = strong["best_overall_delta_E"].values
    left_y = np.arange(len(strong))
    axes[0].barh(left_y, left_vals, color="#9467bd")
    axes[0].set_yticks(left_y)
    axes[0].set_yticklabels(strong["pair"], fontsize=9)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Best delta_E")
    axes[0].set_title(f"Largest-artifact {len(strong)} metal pairs")
    apply_delta_xlim(axes[0], left_vals)

    right_vals = weak["best_overall_delta_E"].values
    right_y = np.arange(len(weak))
    axes[1].barh(right_y, right_vals, color="#8c564b")
    axes[1].set_yticks(right_y)
    axes[1].set_yticklabels(weak["pair"], fontsize=9)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Best delta_E")
    axes[1].set_title(f"Smallest-artifact {len(weak)} metal pairs")
    apply_delta_xlim(axes[1], right_vals)

    save_fig(fig, RANKING_DIR / "pair_extremes_top_vs_bottom.png")


def plot_delta_vs_dip_length(df: pd.DataFrame):
    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.scatter(
        df["dip_cell_length"],
        df["delta_E"],
        s=26,
        color="lightgray",
        alpha=0.65,
        label="All accepted configurations",
    )
    ax.scatter(
        top_cfg["dip_cell_length"],
        top_cfg["delta_E"],
        s=70,
        color="#1f77b4",
        label="Largest artifacts",
    )
    ax.scatter(
        bottom_cfg["dip_cell_length"],
        bottom_cfg["delta_E"],
        s=70,
        color="#d62728",
        label="Smallest accepted artifacts",
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
    ax.set_title("Ternary systems: largest and smallest accepted artifact cases")
    ax.legend()

    vals = pd.concat([df["delta_E"], pd.Series([0.0])], ignore_index=True)
    vals = pd.to_numeric(vals, errors="coerce").dropna()
    ymin = float(vals.min()) if not vals.empty else -1.0
    ymax = float(vals.max()) if not vals.empty else 0.02
    pad = max(0.03, 0.08 * max(abs(ymin), abs(ymax), 1e-6))
    ax.set_ylim(ymin - pad, ymax + pad)

    save_fig(fig, SCATTER_DIR / "delta_vs_dip_cell_length_top_bottom.png")


def plot_delta_vs_ie(df: pd.DataFrame, x_col: str, out_name: str, title: str, xlabel: str):
    top_cfg = top_configurations(df, TOP_N_CONFIGS)
    bottom_cfg = bottom_configurations(df, TOP_N_CONFIGS)

    fig, ax = plt.subplots(figsize=(11, 7))
    ax.scatter(
        df[x_col],
        df["delta_E"],
        s=26,
        color="lightgray",
        alpha=0.65,
        label="All accepted configurations",
    )
    ax.scatter(
        top_cfg[x_col],
        top_cfg["delta_E"],
        s=70,
        color="#1f77b4",
        label="Largest artifacts",
    )
    ax.scatter(
        bottom_cfg[x_col],
        bottom_cfg["delta_E"],
        s=70,
        color="#d62728",
        label="Smallest accepted artifacts",
    )

    for _, row in top_cfg.head(min(ANNOTATE_POINTS, len(top_cfg))).iterrows():
        ax.annotate(
            row["pair"],
            (row[x_col], row["delta_E"]),
            fontsize=8,
            xytext=(4, 3),
            textcoords="offset points",
        )

    for _, row in bottom_cfg.head(min(ANNOTATE_POINTS, len(bottom_cfg))).iterrows():
        ax.annotate(
            row["pair"],
            (row[x_col], row["delta_E"]),
            fontsize=8,
            xytext=(4, -10),
            textcoords="offset points",
        )

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("delta_E")
    ax.set_title(title)
    ax.legend()

    vals = pd.concat([df["delta_E"], pd.Series([0.0])], ignore_index=True)
    vals = pd.to_numeric(vals, errors="coerce").dropna()
    ymin = float(vals.min()) if not vals.empty else -1.0
    ymax = float(vals.max()) if not vals.empty else 0.02
    pad = max(0.03, 0.08 * max(abs(ymin), abs(ymax), 1e-6))
    ax.set_ylim(ymin - pad, ymax + pad)

    save_fig(fig, SCATTER_DIR / out_name)


def plot_ligand_top_bottom_counts(df: pd.DataFrame):
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
    ax.set_title("Ligand makeup of ternary top vs bottom artifact cases")
    ax.legend()

    save_fig(fig, SUMMARY_DIR / "ligand_top_bottom_counts.png")


def plot_pair_heatmap_extremes(df: pd.DataFrame, mode: str):
    summary = pair_summary(df)
    if mode == "top":
        chosen = strongest_pairs(summary, TOP_N_PAIRS)
        value_col = "best_overall_delta_E"
        title = "Top-pair heatmap (largest artifacts)"
        out_name = "top_pair_heatmap.png"
    else:
        chosen = weakest_pairs(summary, TOP_N_PAIRS)
        value_col = "best_overall_delta_E"
        title = "Bottom-pair heatmap (smallest artifacts)"
        out_name = "bottom_pair_heatmap.png"

    pair_set = set(chosen["pair"])
    usable = summary[summary["pair"].isin(pair_set)].copy()

    matrix = pd.DataFrame(np.nan, index=METAL_ORDER, columns=METAL_ORDER)
    grouped = usable.groupby("pair", as_index=False)["best_delta_E"].min()
    for _, row in grouped.iterrows():
        m1, m2 = row["pair"].split("-")
        val = row["best_delta_E"]
        if m1 in matrix.index and m2 in matrix.columns:
            matrix.loc[m1, m2] = val
            matrix.loc[m2, m1] = val

    fig, ax = plt.subplots(figsize=(12, 10))
    arr = matrix.values.astype(float)
    im = ax.imshow(arr, aspect="auto")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Best delta_E")
    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_xticklabels(matrix.columns, rotation=90)
    ax.set_yticklabels(matrix.index)
    ax.set_title(title)

    save_fig(fig, HEATMAP_DIR / out_name)


# ============================================================
# MAIN
# ============================================================
def main():
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input CSV not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)
    if df.empty:
        raise ValueError("Input CSV is empty.")

    required = {
        "chemsys", "metal_1", "metal_2", "ligand",
        "source", "delta_E", "dip_cell_length", "baseline_slope"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = add_descriptor_columns(df)

    numeric_cols = ["delta_E", "dip_cell_length", "baseline_slope", "IE_1", "IE_2", "IE_avg", "IE_diff"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["delta_E", "dip_cell_length"]).copy()
    if df.empty:
        raise ValueError("No rows remain after numeric cleanup.")

    print(f"Loaded rows: {len(df)}")
    print(f"Unique chemsys: {df['chemsys'].nunique()}")
    print(f"Unique pairs: {df['pair'].nunique()}")
    print(f"Ligands: {sorted(df['ligand'].dropna().unique())}")

    write_tables(df)
    plot_configuration_extremes(df)
    plot_chemsys_extremes(df)
    plot_pair_extremes(df)
    plot_delta_vs_dip_length(df)
    plot_delta_vs_ie(
        df,
        "IE_avg",
        "delta_vs_IE_avg_top_bottom.png",
        "Ternary systems: delta_E vs average ionization energy",
        "Average first ionization energy (kJ/mol)",
    )
    plot_delta_vs_ie(
        df,
        "IE_diff",
        "delta_vs_IE_diff_top_bottom.png",
        "Ternary systems: delta_E vs ionization-energy mismatch",
        "|IE_1 - IE_2| (kJ/mol)",
    )
    plot_ligand_top_bottom_counts(df)
    plot_pair_heatmap_extremes(df, "top")
    plot_pair_heatmap_extremes(df, "bottom")

    print("\nDone.")
    print(f"Saved ternary top-vs-bottom analysis to: {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()

