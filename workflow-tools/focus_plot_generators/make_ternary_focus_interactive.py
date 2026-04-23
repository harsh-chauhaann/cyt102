#!/usr/bin/env python3
"""
Create interactive ternary metal-pair focus plots as zoomable HTML files.

Features:
- One HTML plot per ternary metal pair
- Hover tooltips instead of always-on labels
- Parallel batch generation with ProcessPoolExecutor

Examples:
    python make_ternary_focus_interactive.py Ag-Fe
    python make_ternary_focus_interactive.py --all
    python make_ternary_focus_interactive.py --all --workers 8
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path

import pandas as pd
import plotly.express as px


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "packing-fraction-review" / "ternary" / "tables" / "packing_delta_joined.csv"
OUT_DIR = ROOT / "packing-fraction-review" / "ternary" / "interactive_plots"

COLORS = {"O": "#1f77b4", "S": "#d62728"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("pair", nargs="?", default=None, help="Ternary metal pair, e.g. Ag-Fe")
    parser.add_argument("--all", action="store_true", help="Generate plots for all ternary metal pairs")
    parser.add_argument("--workers", type=int, default=None, help="Number of worker processes for --all mode")
    return parser.parse_args()


def canonical_pair(m1: str, m2: str) -> str:
    a, b = sorted([str(m1), str(m2)], key=lambda x: x.upper())
    return f"{a}-{b}"


def pair_list(df: pd.DataFrame) -> list[str]:
    return sorted(df["pair"].dropna().astype(str).unique())


def prepare_df() -> pd.DataFrame:
    df = pd.read_csv(INPUT_CSV)
    if "pair" not in df.columns:
        df["pair"] = df.apply(lambda row: canonical_pair(row["metal_1"], row["metal_2"]), axis=1)
    return df


def prepare_subset(df: pd.DataFrame, pair: str) -> pd.DataFrame:
    subset = df[df["pair"].astype(str) == pair].copy()
    if subset.empty:
        raise ValueError(f"No ternary rows found for pair: {pair}")

    subset["short_source"] = subset["source"].str.replace("_POSCAR", "", regex=False)
    subset = subset.sort_values(["ligand", "packing_fraction", "delta_E"]).reset_index(drop=True)
    return subset


def build_figure(subset: pd.DataFrame, pair: str):
    fig = px.scatter(
        subset,
        x="packing_fraction",
        y="delta_E",
        color="ligand",
        color_discrete_map=COLORS,
        hover_name="short_source",
        hover_data={
            "pair": True,
            "chemsys": True,
            "metal_1": True,
            "metal_2": True,
            "packing_fraction": ":.4f",
            "delta_E": ":.4f",
            "dip_cell_length": ":.4f",
            "formula": True,
            "material_id": True,
            "ligand": False,
            "short_source": False,
        },
        title=f"Ternary {pair} systems: packing efficiency vs delta_E",
    )

    fig.update_traces(marker=dict(size=11, line=dict(width=0.7, color="black")), opacity=0.85)
    fig.update_layout(
        template="plotly_white",
        legend_title_text="Ligand",
        xaxis_title="Packing fraction (covalent-radius proxy)",
        yaxis_title="delta_E",
        hoverlabel=dict(bgcolor="white"),
    )
    fig.add_hline(y=0.0, line_dash="dash", line_color="black")
    fig.update_xaxes(showgrid=True, gridcolor="rgba(0,0,0,0.08)", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(0,0,0,0.08)", zeroline=False)
    return fig


def write_one_plot(task: tuple[str, list[dict]]) -> str:
    pair, records = task
    subset = pd.DataFrame.from_records(records)
    fig = build_figure(subset, pair)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{pair.lower()}_pair_packing_vs_delta_interactive.html"
    fig.write_html(out_path, include_plotlyjs="cdn", full_html=True)
    return str(out_path)


def write_index(pairs: list[str]) -> Path:
    lines = [
        "<html><head><meta charset='utf-8'><title>Ternary Interactive Plots</title></head><body>",
        "<h1>Ternary Interactive Packing-vs-Delta Plots</h1>",
        "<p>Open a metal-pair plot below. Use mouse wheel / drag to zoom and pan. Hover on points to see structure details.</p>",
        "<ul>",
    ]
    for pair in pairs:
        filename = f"{pair.lower()}_pair_packing_vs_delta_interactive.html"
        lines.append(f"<li><a href='{filename}'>{pair}</a></li>")
    lines.extend(["</ul>", "</body></html>"])
    index_path = OUT_DIR / "index.html"
    index_path.write_text("\n".join(lines), encoding="utf-8")
    return index_path


def main() -> None:
    args = parse_args()
    if not args.all and not args.pair:
        raise SystemExit("Pass a metal pair like 'Ag-Fe' or use --all")

    df = prepare_df()
    pairs = pair_list(df)

    if args.all:
        tasks = []
        for pair in pairs:
            subset = prepare_subset(df, pair)
            tasks.append((pair, subset.to_dict(orient="records")))

        max_workers = args.workers or min(len(tasks), max((os.cpu_count() or 2) - 1, 1))
        written: list[str] = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(write_one_plot, task): task[0] for task in tasks}
            for future in as_completed(futures):
                pair = futures[future]
                out_path = future.result()
                written.append(out_path)
                print(f"Saved {pair}: {out_path}")

        index_path = write_index(pairs)
        print(f"Wrote {len(written)} interactive plots.")
        print(f"Index: {index_path}")
        return

    pair = args.pair.strip()
    subset = prepare_subset(df, pair)
    out_path = write_one_plot((pair, subset.to_dict(orient="records")))
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

