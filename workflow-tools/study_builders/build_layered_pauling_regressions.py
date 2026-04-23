#!/usr/bin/env python3
"""
Create a cleaner set of binary electronegativity-vs-delta_E regression plots.

This version is designed for discussion and reporting:
- no per-point text labels
- lighter scatter points
- clearer regression overlays
- a 2x2 statistic comparison grid for O vs S
- an auto-written findings seed markdown file
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
OUT_DIR = ROOT / "layered-pauling-regressions"
PLOTS_DIR = OUT_DIR / "plots"
TABLES_DIR = OUT_DIR / "tables"

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

STAT_CONFIG = {
    "min": {"label": "Min", "color": "#1f77b4", "marker": "o"},
    "max": {"label": "Max", "color": "#d62728", "marker": "s"},
    "mean": {"label": "Mean", "color": "#2ca02c", "marker": "^"},
    "median": {"label": "Median", "color": "#9467bd", "marker": "D"},
}

LIGAND_CONFIG = {
    "O": {"label": "Oxygen", "color": "#1f77b4"},
    "S": {"label": "Sulfur", "color": "#d62728"},
}


def ensure_dirs() -> None:
    for path in [OUT_DIR, PLOTS_DIR, TABLES_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def regression_metrics(x: pd.Series, y: pd.Series) -> dict[str, object] | None:
    valid = x.notna() & y.notna()
    if valid.sum() < 2:
        return None

    x_ok = x[valid].astype(float)
    y_ok = y[valid].astype(float)
    if x_ok.nunique() < 2 or y_ok.nunique() < 2:
        return None

    slope, intercept = np.polyfit(x_ok, y_ok, 1)
    corr = float(x_ok.corr(y_ok))
    x_line = np.linspace(float(x_ok.min()), float(x_ok.max()), 100)
    y_line = slope * x_line + intercept
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "correlation": corr,
        "n_points": int(valid.sum()),
        "x_line": x_line,
        "y_line": y_line,
    }


def prepare_df() -> pd.DataFrame:
    df = pd.read_csv(INPUT_CSV)
    df["electronegativity"] = df["metal"].map(PAULING_ELECTRONEGATIVITY)
    missing = df[df["electronegativity"].isna()]["metal"].drop_duplicates().tolist()
    if missing:
        raise ValueError(f"Missing electronegativity values for: {missing}")
    return df


def overall_y_limits(df: pd.DataFrame) -> tuple[float, float]:
    cols = [f"delta_E_{key}" for key in STAT_CONFIG]
    vals = pd.concat([pd.to_numeric(df[col], errors="coerce") for col in cols], ignore_index=True).dropna()
    ymin = float(vals.min())
    ymax = float(vals.max())
    pad = 0.06 * max(ymax - ymin, 1e-6)
    return ymin - pad, ymax + pad


def write_ligand_clean_plot(
    df: pd.DataFrame,
    ligand: str,
    y_limits: tuple[float, float],
) -> tuple[Path, list[dict[str, object]]]:
    subset = df[df["ligand"] == ligand].copy()
    subset = subset.sort_values(["electronegativity", "metal"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10.2, 6.8))
    summary_rows: list[dict[str, object]] = []

    for stat_key, style in STAT_CONFIG.items():
        stat_col = f"delta_E_{stat_key}"

        ax.scatter(
            subset["electronegativity"],
            subset[stat_col],
            s=52,
            alpha=0.55,
            color=style["color"],
            marker=style["marker"],
            edgecolors="white",
            linewidths=0.4,
            zorder=3,
        )

        fit = regression_metrics(subset["electronegativity"], subset[stat_col])
        if fit is None:
            continue

        ax.plot(
            fit["x_line"],
            fit["y_line"],
            color=style["color"],
            linewidth=2.6,
            linestyle="--",
            label=f"{style['label']} (r={fit['correlation']:.3f})",
            zorder=4,
        )

        summary_rows.append({
            "view": "ligand_overview",
            "ligand": ligand,
            "statistic": stat_key,
            "n_points": fit["n_points"],
            "slope": fit["slope"],
            "intercept": fit["intercept"],
            "correlation": fit["correlation"],
        })

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1.2)
    ax.set_ylim(*y_limits)
    ax.set_xlabel("Metal electronegativity (Pauling)")
    ax.set_ylabel("delta_E")
    ax.set_title(f"Binary {LIGAND_CONFIG[ligand]['label']}: multiple regression overview")
    ax.grid(alpha=0.18)
    ax.legend(title="Statistic", fontsize=9)

    out_path = PLOTS_DIR / f"binary_{ligand.lower()}_clean_multiple_regression.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path, summary_rows


def write_statistic_grid(
    df: pd.DataFrame,
    y_limits: tuple[float, float],
) -> tuple[Path, list[dict[str, object]]]:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    rows: list[dict[str, object]] = []

    for ax, stat_key in zip(axes, STAT_CONFIG):
        stat_col = f"delta_E_{stat_key}"
        for ligand, ligand_style in LIGAND_CONFIG.items():
            subset = df[df["ligand"] == ligand].copy()
            subset = subset.sort_values(["electronegativity", "metal"]).reset_index(drop=True)

            ax.scatter(
                subset["electronegativity"],
                subset[stat_col],
                s=48,
                alpha=0.5,
                color=ligand_style["color"],
                edgecolors="white",
                linewidths=0.35,
                zorder=3,
            )

            fit = regression_metrics(subset["electronegativity"], subset[stat_col])
            if fit is None:
                continue

            ax.plot(
                fit["x_line"],
                fit["y_line"],
                color=ligand_style["color"],
                linewidth=2.5,
                linestyle="-",
                label=f"{ligand} (r={fit['correlation']:.3f})",
                zorder=4,
            )

            rows.append({
                "view": "statistic_grid",
                "ligand": ligand,
                "statistic": stat_key,
                "n_points": fit["n_points"],
                "slope": fit["slope"],
                "intercept": fit["intercept"],
                "correlation": fit["correlation"],
            })

        ax.axhline(0.0, color="black", linestyle="--", linewidth=1.0)
        ax.set_ylim(*y_limits)
        ax.set_title(f"{STAT_CONFIG[stat_key]['label']} delta_E")
        ax.grid(alpha=0.18)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("delta_E")
    axes[2].set_ylabel("delta_E")
    axes[2].set_xlabel("Metal electronegativity (Pauling)")
    axes[3].set_xlabel("Metal electronegativity (Pauling)")
    fig.suptitle("Binary electronegativity trends by statistic", fontsize=15)

    out_path = PLOTS_DIR / "binary_clean_statistic_comparison_grid.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path, rows


def write_findings_seed(summary_df: pd.DataFrame) -> None:
    strongest = (
        summary_df.assign(abs_correlation=summary_df["correlation"].abs())
        .sort_values("abs_correlation", ascending=False)
        .reset_index(drop=True)
    )

    lines = [
        "# Electronegativity Findings Seed",
        "",
        "This file is a starting point for the markdown discussion of the cleaner regression plots.",
        "",
        "## Main observations",
        "",
    ]

    for _, row in strongest.head(4).iterrows():
        direction = "positive" if row["correlation"] > 0 else "negative"
        lines.append(
            f"- `{row['ligand']}` with `{row['statistic']}` shows a `{direction}` correlation "
            f"of `{row['correlation']:.3f}` and slope `{row['slope']:.4f}`."
        )

    lines.extend([
        "",
        "## Suggested discussion points",
        "",
        "- The strongest trends appear in the extreme statistics (`min` or `max`) rather than in `mean` or `median`.",
        "- Oxygen and sulfur do not show the same regression behavior for every statistic.",
        "- Mean and median trends are comparatively weak, so the extreme points may be driving most of the visible behavior.",
        "- Any interpretation should mention the role of outliers, especially for the `min` delta_E values.",
        "",
        "## Files to reference",
        "",
        "- `plots/binary_o_clean_multiple_regression.png`",
        "- `plots/binary_s_clean_multiple_regression.png`",
        "- `plots/binary_clean_statistic_comparison_grid.png`",
        "- `tables/layered_regression_summary.csv`",
    ])

    (OUT_DIR / "drafting_outline.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readme(summary_df: pd.DataFrame) -> None:
    lines = [
        "# Binary Electronegativity Clean Multi-Regression Plots",
        "",
        "Cleaner follow-up plots for discussing electronegativity trends.",
        "",
        "Outputs:",
        "- `plots/binary_o_clean_multiple_regression.png`",
        "- `plots/binary_s_clean_multiple_regression.png`",
        "- `plots/binary_clean_statistic_comparison_grid.png`",
        "- `tables/layered_regression_summary.csv`",
        "- `drafting_outline.md`",
        "",
        f"- regression rows saved: `{len(summary_df)}`",
    ]
    (OUT_DIR / "overview_notes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    df = prepare_df()
    y_limits = overall_y_limits(df)

    all_rows: list[dict[str, object]] = []

    for ligand in ["O", "S"]:
        out_path, rows = write_ligand_clean_plot(df, ligand, y_limits)
        print(f"Saved: {out_path}")
        all_rows.extend(rows)

    grid_path, grid_rows = write_statistic_grid(df, y_limits)
    print(f"Saved: {grid_path}")
    all_rows.extend(grid_rows)

    summary_df = pd.DataFrame(all_rows)
    summary_df["view"] = summary_df["view"].fillna("ligand_overview")
    summary_df = summary_df.sort_values(["view", "ligand", "statistic"]).reset_index(drop=True)
    summary_path = TABLES_DIR / "layered_regression_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved: {summary_path}")

    write_findings_seed(summary_df[summary_df["view"] == "statistic_grid"].copy())
    write_readme(summary_df)


if __name__ == "__main__":
    main()

