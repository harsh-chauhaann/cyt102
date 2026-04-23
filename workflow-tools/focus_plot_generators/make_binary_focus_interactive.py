#!/usr/bin/env python3
"""
Create interactive binary metal focus plots as zoomable HTML files.

Features:
- One HTML plot per binary metal
- Hover tooltips instead of always-on labels
- Parallel batch generation with ProcessPoolExecutor

Examples:
    python make_binary_focus_interactive.py Sc
    python make_binary_focus_interactive.py --all
    python make_binary_focus_interactive.py --all --workers 8
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path

import pandas as pd
import plotly.express as px


ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "packing-fraction-review" / "binary" / "tables" / "packing_delta_joined.csv"
OUT_DIR = ROOT / "packing-fraction-review" / "binary" / "interactive_plots"

COLORS = {"O": "#1f77b4", "S": "#d62728"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("metal", nargs="?", default=None, help="Binary metal symbol, e.g. Sc")
    parser.add_argument("--all", action="store_true", help="Generate plots for all binary metals")
    parser.add_argument("--workers", type=int, default=None, help="Number of worker processes for --all mode")
    return parser.parse_args()


def sorted_metals(df: pd.DataFrame) -> list[str]:
    return sorted({str(value).split("-")[0] for value in df["chemsys"].dropna()})


def prepare_subset(df: pd.DataFrame, metal: str) -> pd.DataFrame:
    subset = df[df["chemsys"].astype(str).str.startswith(f"{metal}-")].copy()
    if subset.empty:
        raise ValueError(f"No binary rows found for metal: {metal}")

    subset["ligand"] = subset["chemsys"].str.split("-").str[1]
    subset["short_source"] = subset["source"].str.replace("_POSCAR", "", regex=False)
    subset = subset.sort_values(["ligand", "packing_fraction", "delta_E"]).reset_index(drop=True)
    return subset


def build_figure(subset: pd.DataFrame, metal: str):
    fig = px.scatter(
        subset,
        x="packing_fraction",
        y="delta_E",
        color="ligand",
        color_discrete_map=COLORS,
        hover_name="short_source",
        hover_data={
            "chemsys": True,
            "packing_fraction": ":.4f",
            "delta_E": ":.4f",
            "dip_cell_length": ":.4f",
            "formula": True,
            "material_id": True,
            "ligand": False,
            "short_source": False,
        },
        title=f"Binary {metal} systems: packing efficiency vs delta_E",
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
    metal, records = task
    subset = pd.DataFrame.from_records(records)
    fig = build_figure(subset, metal)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{metal.lower()}_metal_packing_vs_delta_interactive.html"
    fig.write_html(out_path, include_plotlyjs="cdn", full_html=True)
    return str(out_path)


def write_index(metals: list[str]) -> Path:
    lines = [
        "<html><head><meta charset='utf-8'><title>Binary Interactive Plots</title></head><body>",
        "<h1>Binary Interactive Packing-vs-Delta Plots</h1>",
        "<p>Open a metal plot below. Use mouse wheel / drag to zoom and pan. Hover on points to see structure details.</p>",
        "<ul>",
    ]
    for metal in metals:
        filename = f"{metal.lower()}_metal_packing_vs_delta_interactive.html"
        lines.append(f"<li><a href='{filename}'>{metal}</a></li>")
    lines.extend(["</ul>", "</body></html>"])
    index_path = OUT_DIR / "index.html"
    index_path.write_text("\n".join(lines), encoding="utf-8")
    return index_path


def main() -> None:
    args = parse_args()
    if not args.all and not args.metal:
        raise SystemExit("Pass a metal symbol like 'Sc' or use --all")

    df = pd.read_csv(INPUT_CSV)
    metals = sorted_metals(df)

    if args.all:
        tasks = []
        for metal in metals:
            subset = prepare_subset(df, metal)
            tasks.append((metal, subset.to_dict(orient="records")))

        max_workers = args.workers or min(len(tasks), max((os.cpu_count() or 2) - 1, 1))
        written: list[str] = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(write_one_plot, task): task[0] for task in tasks}
            for future in as_completed(futures):
                metal = futures[future]
                out_path = future.result()
                written.append(out_path)
                print(f"Saved {metal}: {out_path}")

        index_path = write_index(metals)
        print(f"Wrote {len(written)} interactive plots.")
        print(f"Index: {index_path}")
        return

    metal = args.metal.strip()
    subset = prepare_subset(df, metal)
    out_path = write_one_plot((metal, subset.to_dict(orient="records")))
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

