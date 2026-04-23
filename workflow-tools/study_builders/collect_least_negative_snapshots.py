#!/usr/bin/env python3
"""
Build summary tables and plots using the maximum delta_E point per category:
- Binary: per metal and ligand (O/S)
- Ternary: per metal pair and ligand (O/S)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "least-negative-delta-snapshots"
TABLES_DIR = OUT_DIR / "tables"
PLOTS_DIR = OUT_DIR / "plots"

BINARY_INPUT = ROOT / "packing-fraction-review" / "binary" / "tables" / "packing_delta_joined.csv"
TERNARY_INPUT = ROOT / "packing-fraction-review" / "ternary" / "tables" / "packing_delta_joined.csv"

COLORS = {"O": "#1f77b4", "S": "#d62728"}


def ensure_dirs() -> None:
    for path in [OUT_DIR, TABLES_DIR, PLOTS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def canonical_pair(m1: str, m2: str) -> str:
    a, b = sorted([str(m1), str(m2)], key=lambda x: x.upper())
    return f"{a}-{b}"


def highest_delta_rows_binary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["metal"] = out["chemsys"].str.split("-").str[0]
    out["ligand"] = out["chemsys"].str.split("-").str[1]
    idx = out.groupby(["metal", "ligand"])["delta_E"].idxmax()
    chosen = out.loc[idx].copy()
    chosen = chosen.sort_values(["metal", "ligand"]).reset_index(drop=True)
    return chosen


def highest_delta_rows_ternary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "pair" not in out.columns:
        out["pair"] = out.apply(lambda row: canonical_pair(row["metal_1"], row["metal_2"]), axis=1)
    idx = out.groupby(["pair", "ligand"])["delta_E"].idxmax()
    chosen = out.loc[idx].copy()
    chosen = chosen.sort_values(["pair", "ligand"]).reset_index(drop=True)
    return chosen


def save_binary_png(df: pd.DataFrame, out_path: Path) -> None:
    metals = sorted(df["metal"].dropna().astype(str).unique())
    pos = np.arange(len(metals))
    offsets = {"O": -0.12, "S": 0.12}

    fig, ax = plt.subplots(figsize=(14, 6))
    for ligand in ["O", "S"]:
        group = df[df["ligand"] == ligand].copy()
        x_vals = [pos[metals.index(metal)] + offsets[ligand] for metal in group["metal"]]
        ax.scatter(
            x_vals,
            group["delta_E"],
            s=85,
            color=COLORS[ligand],
            edgecolors="black",
            linewidths=0.6,
            label=ligand,
            zorder=3,
        )

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(pos)
    ax.set_xticklabels(metals)
    ax.set_ylabel("Maximum delta_E")
    ax.set_xlabel("Metal")
    ax.set_title("Binary: maximum delta_E point per metal for O and S")
    ax.grid(axis="y", alpha=0.18)
    ax.legend(title="Ligand")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_ternary_png(df: pd.DataFrame, out_path: Path) -> None:
    pairs = sorted(df["pair"].dropna().astype(str).unique())
    pos = np.arange(len(pairs))
    offsets = {"O": -0.12, "S": 0.12}

    fig_w = max(24, 0.22 * len(pairs))
    fig, ax = plt.subplots(figsize=(fig_w, 7))
    for ligand in ["O", "S"]:
        group = df[df["ligand"] == ligand].copy()
        x_vals = [pos[pairs.index(pair)] + offsets[ligand] for pair in group["pair"]]
        ax.scatter(
            x_vals,
            group["delta_E"],
            s=55,
            color=COLORS[ligand],
            edgecolors="black",
            linewidths=0.5,
            label=ligand,
            zorder=3,
        )

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(pos)
    ax.set_xticklabels(pairs, rotation=90, fontsize=7)
    ax.set_ylabel("Maximum delta_E")
    ax.set_xlabel("Metal pair")
    ax.set_title("Ternary: maximum delta_E point per metal pair for O and S")
    ax.grid(axis="y", alpha=0.18)
    ax.legend(title="Ligand")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_binary_html(df: pd.DataFrame, out_path: Path) -> None:
    fig = px.scatter(
        df,
        x="metal",
        y="delta_E",
        color="ligand",
        color_discrete_map=COLORS,
        hover_name="source",
        hover_data={
            "chemsys": True,
            "delta_E": ":.4f",
            "dip_cell_length": ":.4f",
            "formula": True,
            "material_id": True,
            "ligand": False,
            "metal": False,
        },
        title="Binary: maximum delta_E point per metal for O and S",
    )
    fig.update_traces(marker=dict(size=11, line=dict(width=0.7, color="black")), opacity=0.9)
    fig.update_layout(template="plotly_white", xaxis_title="Metal", yaxis_title="Maximum delta_E")
    fig.add_hline(y=0.0, line_dash="dash", line_color="black")
    fig.write_html(out_path, include_plotlyjs="cdn", full_html=True)


def save_ternary_html(df: pd.DataFrame, out_path: Path) -> None:
    fig = px.scatter(
        df,
        x="pair",
        y="delta_E",
        color="ligand",
        color_discrete_map=COLORS,
        hover_name="source",
        hover_data={
            "pair": True,
            "chemsys": True,
            "delta_E": ":.4f",
            "dip_cell_length": ":.4f",
            "formula": True,
            "material_id": True,
            "ligand": False,
        },
        title="Ternary: maximum delta_E point per metal pair for O and S",
    )
    fig.update_traces(marker=dict(size=9, line=dict(width=0.6, color="black")), opacity=0.9)
    fig.update_layout(
        template="plotly_white",
        xaxis_title="Metal pair",
        yaxis_title="Maximum delta_E",
        xaxis=dict(tickangle=90),
    )
    fig.add_hline(y=0.0, line_dash="dash", line_color="black")
    fig.write_html(out_path, include_plotlyjs="cdn", full_html=True)


def write_readme(binary_df: pd.DataFrame, ternary_df: pd.DataFrame) -> None:
    lines = [
        "# Max delta_E Point Results",
        "",
        "These outputs keep only the single highest `delta_E` row in each category.",
        "",
        "Interpretation used here:",
        "- `maximum delta_E` means the numerically largest `delta_E` value",
        "- for binary, categories are `(metal, ligand)`",
        "- for ternary, categories are `(metal pair, ligand)`",
        "",
        f"- binary selected rows: `{len(binary_df)}`",
        f"- ternary selected rows: `{len(ternary_df)}`",
        "",
        "Files:",
        "- `tables/binary_least_negative_points.csv`",
        "- `tables/ternary_least_negative_points.csv`",
        "- `plots/binary_max_delta_points.png`",
        "- `plots/ternary_max_delta_points.png`",
        "- `plots/binary_max_delta_points.html`",
        "- `plots/ternary_max_delta_points.html`",
    ]
    (OUT_DIR / "overview_notes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    binary_df = pd.read_csv(BINARY_INPUT)
    ternary_df = pd.read_csv(TERNARY_INPUT)

    binary_max = highest_delta_rows_binary(binary_df)
    ternary_max = highest_delta_rows_ternary(ternary_df)

    binary_max.to_csv(TABLES_DIR / "binary_least_negative_points.csv", index=False)
    ternary_max.to_csv(TABLES_DIR / "ternary_least_negative_points.csv", index=False)

    save_binary_png(binary_max, PLOTS_DIR / "binary_max_delta_points.png")
    save_ternary_png(ternary_max, PLOTS_DIR / "ternary_max_delta_points.png")
    save_binary_html(binary_max, PLOTS_DIR / "binary_max_delta_points.html")
    save_ternary_html(ternary_max, PLOTS_DIR / "ternary_max_delta_points.html")
    write_readme(binary_max, ternary_max)

    print(f"Saved results to: {OUT_DIR}")


if __name__ == "__main__":
    main()

