#!/usr/bin/env python3
"""
Create bar graphs for grouped delta_E statistics:
- Binary O by metal
- Binary S by metal
- Ternary O by metal pair
- Ternary S by metal pair

For each of the above, generate min / max / mean / median versions and print the
numeric value near the bar tip.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BINARY_INPUT = ROOT / "packing-fraction-review" / "binary" / "tables" / "packing_delta_joined.csv"
TERNARY_INPUT = ROOT / "packing-fraction-review" / "ternary" / "tables" / "packing_delta_joined.csv"

OUT_DIR = ROOT / "aggregate-comparison-bars"
TABLES_DIR = OUT_DIR / "tables"
PLOTS_DIR = OUT_DIR / "plots"

COLORS = {"O": "#1f77b4", "S": "#d62728"}
STAT_FUNCS = {"min": "min", "max": "max", "mean": "mean", "median": "median"}


def ensure_dirs() -> None:
    for path in [OUT_DIR, TABLES_DIR, PLOTS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def canonical_pair(m1: str, m2: str) -> str:
    a, b = sorted([str(m1), str(m2)], key=lambda x: x.upper())
    return f"{a}-{b}"


def summarize_binary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["metal"] = out["chemsys"].str.split("-").str[0]
    out["ligand"] = out["chemsys"].str.split("-").str[1]
    summary = (
        out.groupby(["metal", "ligand"], as_index=False)
        .agg(
            n_configurations=("delta_E", "size"),
            delta_E_min=("delta_E", "min"),
            delta_E_max=("delta_E", "max"),
            delta_E_mean=("delta_E", "mean"),
            delta_E_median=("delta_E", "median"),
        )
        .sort_values(["metal", "ligand"])
        .reset_index(drop=True)
    )
    return summary


def summarize_ternary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "pair" not in out.columns:
        out["pair"] = out.apply(lambda row: canonical_pair(row["metal_1"], row["metal_2"]), axis=1)
    summary = (
        out.groupby(["pair", "ligand"], as_index=False)
        .agg(
            n_configurations=("delta_E", "size"),
            delta_E_min=("delta_E", "min"),
            delta_E_max=("delta_E", "max"),
            delta_E_mean=("delta_E", "mean"),
            delta_E_median=("delta_E", "median"),
        )
        .sort_values(["pair", "ligand"])
        .reset_index(drop=True)
    )
    return summary


def add_value_labels(ax, bars, values: pd.Series, rotation: int = 0, fontsize: float = 8.0) -> None:
    valid = pd.to_numeric(values, errors="coerce").dropna()
    span = max(float(valid.max()) - float(valid.min()), 1e-6) if not valid.empty else 1.0
    pad = 0.03 * span
    if pad < 0.01:
        pad = 0.01

    for bar, value in zip(bars, values):
        x = bar.get_x() + bar.get_width() / 2.0
        y = float(value)
        label = f"{y:.3f}"
        if y <= 0:
            ax.text(
                x,
                y - pad,
                label,
                ha="center",
                va="top",
                rotation=rotation,
                fontsize=fontsize,
                color="#222222",
            )
        else:
            ax.text(
                x,
                y + pad,
                label,
                ha="center",
                va="bottom",
                rotation=rotation,
                fontsize=fontsize,
                color="#222222",
            )


def set_y_limits(ax, values: pd.Series) -> None:
    valid = pd.to_numeric(values, errors="coerce").dropna()
    if valid.empty:
        return
    ymin = float(valid.min())
    ymax = float(valid.max())
    span = max(ymax - ymin, 1e-6)
    bottom_pad = max(0.08, 0.16 * span)
    top_pad = max(0.03, 0.08 * span)
    ax.set_ylim(ymin - bottom_pad, ymax + top_pad)


def save_bar_plot(
    df: pd.DataFrame,
    category_col: str,
    ligand: str,
    stat_key: str,
    title: str,
    xlabel: str,
    out_path: Path,
    rotate_ticks: int = 0,
    width_scale: float = 0.35,
    min_width: float = 10.0,
    height: float = 6.0,
    value_rotation: int = 0,
    value_fontsize: float = 8.0,
) -> None:
    stat_col = f"delta_E_{stat_key}"
    subset = df[df["ligand"] == ligand].copy()
    subset = subset.sort_values(category_col).reset_index(drop=True)

    fig_w = max(min_width, width_scale * max(len(subset), 1))
    fig, ax = plt.subplots(figsize=(fig_w, height))
    bars = ax.bar(
        subset[category_col],
        subset[stat_col],
        color=COLORS[ligand],
        edgecolor="black",
        linewidth=0.6,
    )

    add_value_labels(ax, bars, subset[stat_col], rotation=value_rotation, fontsize=value_fontsize)
    set_y_limits(ax, subset[stat_col])

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"{stat_key.capitalize()} delta_E")
    ax.grid(axis="y", alpha=0.18)
    if rotate_ticks:
        ax.tick_params(axis="x", rotation=rotate_ticks, labelsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_readme(binary_summary: pd.DataFrame, ternary_summary: pd.DataFrame) -> None:
    lines = [
        "# Bar Graph Statistics",
        "",
        "These plots summarize `delta_E` by grouped categories.",
        "",
        "Binary grouping:",
        "- category = `(metal, ligand)`",
        "",
        "Ternary grouping:",
        "- category = `(metal pair, ligand)`",
        "",
        "Statistics shown:",
        "- `min`",
        "- `max`",
        "- `mean`",
        "- `median`",
        "",
        "Each plot prints the numeric value near the bar tip.",
        "",
        f"- binary grouped rows: `{len(binary_summary)}`",
        f"- ternary grouped rows: `{len(ternary_summary)}`",
    ]
    (OUT_DIR / "overview_notes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    binary_df = pd.read_csv(BINARY_INPUT)
    ternary_df = pd.read_csv(TERNARY_INPUT)

    binary_summary = summarize_binary(binary_df)
    ternary_summary = summarize_ternary(ternary_df)

    binary_summary.to_csv(TABLES_DIR / "binary_stat_summary.csv", index=False)
    ternary_summary.to_csv(TABLES_DIR / "ternary_stat_summary.csv", index=False)

    for stat_key in STAT_FUNCS:
        save_bar_plot(
            binary_summary,
            category_col="metal",
            ligand="O",
            stat_key=stat_key,
            title=f"Binary metals: {stat_key} delta_E for oxygen systems",
            xlabel="Metal",
            out_path=PLOTS_DIR / f"binary_oxygen_{stat_key}_delta_bar.png",
        )
        save_bar_plot(
            binary_summary,
            category_col="metal",
            ligand="S",
            stat_key=stat_key,
            title=f"Binary metals: {stat_key} delta_E for sulfur systems",
            xlabel="Metal",
            out_path=PLOTS_DIR / f"binary_sulfur_{stat_key}_delta_bar.png",
        )
        save_bar_plot(
            ternary_summary,
            category_col="pair",
            ligand="O",
            stat_key=stat_key,
            title=f"Ternary metal pairs: {stat_key} delta_E for oxygen systems",
            xlabel="Metal pair",
            out_path=PLOTS_DIR / f"ternary_oxygen_{stat_key}_delta_bar.png",
            rotate_ticks=90,
            width_scale=0.28,
            min_width=28.0,
            height=7.5,
            value_rotation=90,
            value_fontsize=5.5,
        )
        save_bar_plot(
            ternary_summary,
            category_col="pair",
            ligand="S",
            stat_key=stat_key,
            title=f"Ternary metal pairs: {stat_key} delta_E for sulfur systems",
            xlabel="Metal pair",
            out_path=PLOTS_DIR / f"ternary_sulfur_{stat_key}_delta_bar.png",
            rotate_ticks=90,
            width_scale=0.28,
            min_width=28.0,
            height=7.5,
            value_rotation=90,
            value_fontsize=5.5,
        )

    write_readme(binary_summary, ternary_summary)
    print(f"Saved grouped bar-graph statistics to: {OUT_DIR}")


if __name__ == "__main__":
    main()

