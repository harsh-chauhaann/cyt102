#!/usr/bin/env python3
"""
Create binary electronegativity-vs-delta_E plots using metal Pauling electronegativity.

Outputs:
- 8 plots total:
  - oxygen: min / max / mean / median
  - sulfur: min / max / mean / median
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "aggregate-comparison-bars" / "tables" / "binary_stat_summary.csv"
OUT_DIR = ROOT / "pauling-trend-checks"
PLOTS_DIR = OUT_DIR / "plots"
TABLES_DIR = OUT_DIR / "tables"

COLORS = {"O": "#1f77b4", "S": "#d62728"}

PAULING_ELECTRONEGATIVITY = {
    "Sc": 1.36,
    "Ti": 1.54,
    "V": 1.63,
    "Cr": 1.66,
    "Mn": 1.55,
    "Fe": 1.83,
    "Co": 1.88,
    "Ni": 1.91,
    "Cu": 1.90,
    "Zn": 1.65,
    "Y": 1.22,
    "Zr": 1.33,
    "Nb": 1.60,
    "Mo": 2.16,
    "Tc": 1.90,
    "Rh": 2.28,
    "Pd": 2.20,
    "Hf": 1.30,
    "Ta": 1.50,
}

STAT_LABELS = {
    "min": "Most negative delta_E",
    "max": "Maximum delta_E",
    "mean": "Mean delta_E",
    "median": "Median delta_E",
}


def ensure_dirs() -> None:
    for path in [OUT_DIR, PLOTS_DIR, TABLES_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def correlation_or_nan(x: pd.Series, y: pd.Series) -> float:
    valid = x.notna() & y.notna()
    if valid.sum() < 2:
        return float("nan")
    x_ok = x[valid]
    y_ok = y[valid]
    if x_ok.nunique() < 2 or y_ok.nunique() < 2:
        return float("nan")
    return float(x_ok.corr(y_ok))


def line_fit(x: pd.Series, y: pd.Series) -> tuple[np.ndarray, np.ndarray] | None:
    valid = x.notna() & y.notna()
    if valid.sum() < 2:
        return None
    x_ok = x[valid].astype(float)
    y_ok = y[valid].astype(float)
    if x_ok.nunique() < 2:
        return None
    slope, intercept = np.polyfit(x_ok, y_ok, 1)
    x_line = np.linspace(float(x_ok.min()), float(x_ok.max()), 100)
    y_line = slope * x_line + intercept
    return x_line, y_line


def plot_one(df: pd.DataFrame, ligand: str, stat_key: str, out_path: Path) -> None:
    stat_col = f"delta_E_{stat_key}"
    subset = df[df["ligand"] == ligand].copy()
    subset = subset.sort_values(["electronegativity", "metal"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    ax.scatter(
        subset["electronegativity"],
        subset[stat_col],
        s=70,
        color=COLORS[ligand],
        edgecolors="black",
        linewidths=0.5,
        alpha=0.9,
        zorder=3,
    )

    fit = line_fit(subset["electronegativity"], subset[stat_col])
    if fit is not None:
        x_line, y_line = fit
        ax.plot(x_line, y_line, color="#444444", linewidth=1.4, linestyle="--", zorder=2)

    offsets = [(5, 5), (5, -12), (-18, 7), (-18, -12), (7, 14), (7, -20)]
    for i, (_, row) in enumerate(subset.iterrows()):
        dx, dy = offsets[i % len(offsets)]
        ax.annotate(
            row["metal"],
            (row["electronegativity"], row[stat_col]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=COLORS[ligand], alpha=0.9),
        )

    corr = correlation_or_nan(subset["electronegativity"], subset[stat_col])
    corr_text = f"corr = {corr:.3f}" if pd.notna(corr) else "corr = NA"

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Metal electronegativity (Pauling)")
    ax.set_ylabel(STAT_LABELS[stat_key])
    ax.set_title(f"Binary {ligand}: {STAT_LABELS[stat_key]} vs metal electronegativity\n{corr_text}")
    ax.grid(alpha=0.18)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_readme(df: pd.DataFrame) -> None:
    metals = ", ".join(sorted(df["metal"].dropna().astype(str).unique()))
    lines = [
        "# Binary Electronegativity vs delta_E",
        "",
        "These plots compare metal electronegativity against grouped binary delta_E statistics.",
        "",
        "Ligands plotted separately:",
        "- O",
        "- S",
        "",
        "Statistics plotted separately:",
        "- min",
        "- max",
        "- mean",
        "- median",
        "",
        f"- metals included: {metals}",
        f"- grouped rows: `{len(df)}`",
    ]
    (OUT_DIR / "overview_notes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    df = pd.read_csv(INPUT_CSV)
    df["electronegativity"] = df["metal"].map(PAULING_ELECTRONEGATIVITY)
    missing = df[df["electronegativity"].isna()]["metal"].drop_duplicates().tolist()
    if missing:
        raise ValueError(f"Missing electronegativity values for: {missing}")

    df.to_csv(TABLES_DIR / "binary_pauling_summary.csv", index=False)

    for ligand, ligand_name in [("O", "oxygen"), ("S", "sulfur")]:
        for stat_key in ["min", "max", "mean", "median"]:
            plot_one(
                df,
                ligand=ligand,
                stat_key=stat_key,
                out_path=PLOTS_DIR / f"binary_{ligand_name}_{stat_key}_electronegativity_vs_delta.png",
            )

    write_readme(df)
    print(f"Saved binary electronegativity plots to: {OUT_DIR}")


if __name__ == "__main__":
    main()

