#!/usr/bin/env python3
"""
Compute packing-efficiency descriptors for configurations that already have delta_E rows,
then save merged tables and comparison plots for binary and ternary datasets.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymatgen.core import Structure


ROOT = Path(__file__).resolve().parents[2]
OUT_ROOT = ROOT / "packing-fraction-review"


COVALENT_RADII = {
    "O": 0.66,
    "S": 1.05,
    "Sc": 1.70,
    "Ti": 1.60,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.39,
    "Fe": 1.32,
    "Co": 1.26,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Y": 1.90,
    "Zr": 1.75,
    "Nb": 1.64,
    "Mo": 1.54,
    "Tc": 1.47,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "Ag": 1.45,
    "Cd": 1.44,
    "Hf": 1.75,
    "Ta": 1.70,
    "W": 1.62,
    "Re": 1.51,
    "Os": 1.44,
    "Ir": 1.41,
    "Pt": 1.36,
    "Au": 1.36,
    "Hg": 1.32,
}


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    raw_dir: Path
    delta_path: Path
    output_dir: Path
    required_delta_cols: tuple[str, ...]


def sanitize_filename(value: object) -> str:
    if value is None:
        return "unknown"
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(value).strip())


def sphere_volume(radius: float) -> float:
    return (4.0 / 3.0) * math.pi * (radius ** 3)


def site_symbol(specie) -> str:
    element = getattr(specie, "element", None)
    if element is not None:
        return element.symbol
    return getattr(specie, "symbol", str(specie))


def packing_metrics_from_cif(cif_text: str) -> dict[str, float]:
    structure = Structure.from_str(cif_text, fmt="cif")
    cell_volume = float(structure.volume)
    n_sites = int(len(structure))

    total_sphere_volume = 0.0
    unknown_species: list[str] = []

    for site in structure:
        for specie, occu in site.species.items():
            symbol = site_symbol(specie)
            radius = COVALENT_RADII.get(symbol)
            if radius is None:
                unknown_species.append(symbol)
                continue
            total_sphere_volume += float(occu) * sphere_volume(radius)

    if unknown_species:
        unknown = ", ".join(sorted(set(unknown_species)))
        raise ValueError(f"Missing covalent radius for species: {unknown}")

    packing_fraction = total_sphere_volume / cell_volume if cell_volume > 0 else np.nan
    return {
        "n_sites": n_sites,
        "cell_volume": cell_volume,
        "volume_per_atom": (cell_volume / n_sites) if n_sites else np.nan,
        "sphere_volume_sum": total_sphere_volume,
        "packing_fraction": packing_fraction,
        "packing_percent": 100.0 * packing_fraction if pd.notna(packing_fraction) else np.nan,
    }


def source_name(formula: object, material_id: object) -> str:
    formula_part = sanitize_filename(formula)
    material_part = sanitize_filename(str(material_id).replace("mp-", "").strip())
    return f"{formula_part}_{material_part}_POSCAR"


def build_target_map(delta_df: pd.DataFrame) -> dict[str, set[str]]:
    targets: dict[str, set[str]] = {}
    for chemsys, group in delta_df.groupby("chemsys"):
        targets[str(chemsys)] = set(group["source"].astype(str))
    return targets


def iter_target_rows(raw_df: pd.DataFrame, wanted_sources: set[str]) -> Iterable[tuple[pd.Series, str]]:
    for _, row in raw_df.iterrows():
        src = source_name(row.get("Formula"), row.get("Material ID"))
        if src in wanted_sources:
            yield row, src


def compute_configuration_metrics(spec: DatasetSpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    delta_df = pd.read_csv(spec.delta_path)
    missing = set(spec.required_delta_cols) - set(delta_df.columns)
    if missing:
        raise ValueError(f"{spec.name}: missing delta columns {sorted(missing)}")

    target_map = build_target_map(delta_df)
    records: list[dict[str, object]] = []
    errors: list[dict[str, object]] = []

    for csv_path in sorted(spec.raw_dir.glob("*.csv")):
        chemsys = csv_path.name.replace("_table_export_with_structures.csv", "")
        wanted_sources = target_map.get(chemsys)
        if not wanted_sources:
            continue

        raw_df = pd.read_csv(csv_path)
        for row, src in iter_target_rows(raw_df, wanted_sources):
            cif_text = row.get("Structure")
            if pd.isna(cif_text) or not str(cif_text).strip():
                errors.append({
                    "chemsys": chemsys,
                    "source": src,
                    "reason": "empty Structure field",
                })
                continue

            try:
                metrics = packing_metrics_from_cif(str(cif_text))
                record = {
                    "chemsys": chemsys,
                    "source": src,
                    "material_id": row.get("Material ID"),
                    "formula": row.get("Formula"),
                    "crystal_system": row.get("Crystal System"),
                    "space_group_symbol": row.get("Space Group Symbol"),
                    "space_group_number": row.get("Space Group Number"),
                    "energy_above_hull": row.get("Energy Above Hull"),
                    "formation_energy": row.get("Formation Energy"),
                    "predicted_stable": row.get("Predicted Stable"),
                    "band_gap": row.get("Band Gap"),
                    "is_metal": row.get("Is Metal"),
                    "density": row.get("Density"),
                }
                record.update(metrics)
                records.append(record)
            except Exception as exc:
                errors.append({
                    "chemsys": chemsys,
                    "source": src,
                    "reason": str(exc),
                })

    metrics_df = pd.DataFrame(records)
    errors_df = pd.DataFrame(errors, columns=["chemsys", "source", "reason"])
    return metrics_df, errors_df


def correlation_or_nan(x: pd.Series, y: pd.Series) -> float:
    valid = x.notna() & y.notna()
    if valid.sum() < 2:
        return np.nan
    x_ok = x[valid]
    y_ok = y[valid]
    if x_ok.nunique() < 2 or y_ok.nunique() < 2:
        return np.nan
    return float(x_ok.corr(y_ok))


def slope_or_nan(x: pd.Series, y: pd.Series) -> float:
    valid = x.notna() & y.notna()
    if valid.sum() < 2:
        return np.nan
    x_ok = x[valid]
    y_ok = y[valid]
    if x_ok.nunique() < 2:
        return np.nan
    return float(np.polyfit(x_ok, y_ok, 1)[0])


def chemsys_summary(merged: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for chemsys, group in merged.groupby("chemsys"):
        rows.append({
            "chemsys": chemsys,
            "n_configurations": len(group),
            "ligand_examples": ",".join(sorted(set(group["ligand"].dropna().astype(str)))) if "ligand" in group.columns else "",
            "packing_fraction_min": group["packing_fraction"].min(),
            "packing_fraction_median": group["packing_fraction"].median(),
            "packing_fraction_max": group["packing_fraction"].max(),
            "packing_fraction_range": group["packing_fraction"].max() - group["packing_fraction"].min(),
            "volume_per_atom_median": group["volume_per_atom"].median(),
            "delta_E_min": group["delta_E"].min(),
            "delta_E_mean": group["delta_E"].mean(),
            "delta_E_max": group["delta_E"].max(),
            "packing_vs_delta_corr": correlation_or_nan(group["packing_fraction"], group["delta_E"]),
            "packing_vs_delta_slope": slope_or_nan(group["packing_fraction"], group["delta_E"]),
        })
    return pd.DataFrame(rows).sort_values(["delta_E_min", "chemsys"], ascending=[True, True]).reset_index(drop=True)


def save_fig(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_configuration_scatter(merged: pd.DataFrame, dataset_name: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    ligands = sorted(merged["ligand"].dropna().astype(str).unique()) if "ligand" in merged.columns else []
    palette = {"O": "#1f77b4", "S": "#d62728"}

    if ligands:
        for ligand in ligands:
            subset = merged[merged["ligand"].astype(str) == ligand]
            ax.scatter(
                subset["packing_fraction"],
                subset["delta_E"],
                s=28,
                alpha=0.7,
                label=ligand,
                color=palette.get(ligand, None),
            )
    else:
        ax.scatter(merged["packing_fraction"], merged["delta_E"], s=28, alpha=0.7, color="#1f77b4")

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Packing fraction (covalent-radius proxy)")
    ax.set_ylabel("delta_E")
    ax.set_title(f"{dataset_name}: configuration-level packing efficiency vs delta_E")
    if ligands:
        ax.legend(title="Ligand")
    save_fig(fig, out_path)


def plot_chemsys_scatter(summary: pd.DataFrame, dataset_name: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(
        summary["packing_fraction_median"],
        summary["delta_E_min"],
        s=38,
        alpha=0.75,
        color="#2ca02c",
    )

    for _, row in summary.head(min(10, len(summary))).iterrows():
        ax.annotate(
            row["chemsys"],
            (row["packing_fraction_median"], row["delta_E_min"]),
            fontsize=8,
            xytext=(4, 3),
            textcoords="offset points",
        )

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Median packing fraction per chemsys")
    ax.set_ylabel("Best (most negative) delta_E in chemsys")
    ax.set_title(f"{dataset_name}: inter-compound comparison")
    save_fig(fig, out_path)


def plot_correlation_hist(summary: pd.DataFrame, dataset_name: str, out_path: Path) -> None:
    usable = summary["packing_vs_delta_corr"].dropna()
    fig, ax = plt.subplots(figsize=(8, 5))
    if usable.empty:
        ax.text(0.5, 0.5, "No chemsys with enough points", ha="center", va="center", transform=ax.transAxes)
    else:
        ax.hist(usable, bins=20, color="#9467bd", edgecolor="black", alpha=0.8)
        ax.axvline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Within-chemsys corr(packing_fraction, delta_E)")
    ax.set_ylabel("Number of chemsys")
    ax.set_title(f"{dataset_name}: intra-compound trend distribution")
    save_fig(fig, out_path)


def write_summary_markdown(
    spec: DatasetSpec,
    merged: pd.DataFrame,
    summary: pd.DataFrame,
    errors: pd.DataFrame,
    out_path: Path,
) -> None:
    overall_corr = correlation_or_nan(merged["packing_fraction"], merged["delta_E"])
    with_corr = summary["packing_vs_delta_corr"].dropna()
    strongest_negative = summary.sort_values("packing_vs_delta_corr", ascending=True).head(10)
    strongest_positive = summary.sort_values("packing_vs_delta_corr", ascending=False).head(10)

    lines = [
        f"# {spec.name.capitalize()} Packing Efficiency vs delta_E",
        "",
        f"- rows in merged table: `{len(merged)}`",
        f"- unique chemsys: `{merged['chemsys'].nunique()}`",
        f"- unique configurations: `{merged['source'].nunique()}`",
        f"- overall corr(packing_fraction, delta_E): `{overall_corr:.4f}`" if pd.notna(overall_corr) else "- overall corr(packing_fraction, delta_E): `NA`",
        f"- chemsys with usable within-chemsys correlation: `{len(with_corr)}`",
        f"- parse / matching errors: `{len(errors)}`",
        "",
        "## Tables",
        "",
        "- `structure_packing_metrics.csv`: one row per configuration with packing descriptors",
        "- `packing_delta_joined.csv`: delta rows merged with packing descriptors",
        "- `system_level_overview.csv`: one row per chemsys for inter-compound comparison",
        "- `packing_match_issues.csv`: parse or matching failures",
        "",
        "## Recommended starting points",
        "",
        "- Use `packing_delta_joined.csv` for configuration-level analysis",
        "- Use `system_level_overview.csv` for equal-weight chemsys comparison",
        "- Use `packing_vs_delta_corr` and `packing_vs_delta_slope` for intra-compound trends",
        "",
        "## Strongest negative within-chemsys correlations",
        "",
    ]

    for _, row in strongest_negative.iterrows():
        value = row["packing_vs_delta_corr"]
        if pd.isna(value):
            continue
        lines.append(f"- `{row['chemsys']}`: `{value:.4f}` with `{int(row['n_configurations'])}` configurations")

    lines.extend([
        "",
        "## Strongest positive within-chemsys correlations",
        "",
    ])

    for _, row in strongest_positive.iterrows():
        value = row["packing_vs_delta_corr"]
        if pd.isna(value):
            continue
        lines.append(f"- `{row['chemsys']}`: `{value:.4f}` with `{int(row['n_configurations'])}` configurations")

    if len(errors):
        lines.extend([
            "",
            "## Notes",
            "",
            "- Some configurations could not be parsed or matched; see `packing_match_issues.csv`.",
        ])

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def ensure_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def analyze_dataset(spec: DatasetSpec) -> None:
    tables_dir = spec.output_dir / "tables"
    plots_dir = spec.output_dir / "plots"
    ensure_dirs([spec.output_dir, tables_dir, plots_dir])

    delta_df = pd.read_csv(spec.delta_path)
    metrics_df, errors_df = compute_configuration_metrics(spec)

    merged = delta_df.merge(metrics_df, on=["chemsys", "source"], how="left", validate="many_to_one")
    merged = merged.dropna(subset=["packing_fraction", "delta_E"]).copy()

    summary_df = chemsys_summary(merged)

    metrics_df.to_csv(tables_dir / "structure_packing_metrics.csv", index=False)
    merged.to_csv(tables_dir / "packing_delta_joined.csv", index=False)
    summary_df.to_csv(tables_dir / "system_level_overview.csv", index=False)
    errors_df.to_csv(tables_dir / "packing_match_issues.csv", index=False)

    plot_configuration_scatter(
        merged,
        spec.name.capitalize(),
        plots_dir / "configuration_level_packing_vs_delta.png",
    )
    plot_chemsys_scatter(
        summary_df,
        spec.name.capitalize(),
        plots_dir / "chemsys_level_median_packing_vs_best_delta.png",
    )
    plot_correlation_hist(
        summary_df,
        spec.name.capitalize(),
        plots_dir / "within_chemsys_correlation_histogram.png",
    )
    write_summary_markdown(
        spec,
        merged,
        summary_df,
        errors_df,
        spec.output_dir / ("binary_notes.md" if spec.name == "binary" else "ternary_notes.md"),
    )


def main() -> None:
    ensure_dirs([OUT_ROOT])
    specs = [
        DatasetSpec(
            name="binary",
            raw_dir=ROOT / "binary-structure-sweep" / "csv",
            delta_path=ROOT / "binary-structure-sweep" / "delta_analysis" / "all_deltaE_combined.csv",
            output_dir=OUT_ROOT / "binary",
            required_delta_cols=("chemsys", "source", "delta_E"),
        ),
        DatasetSpec(
            name="ternary",
            raw_dir=ROOT / "ternary-structure-sweep" / "csv",
            delta_path=ROOT / "ternary-structure-sweep" / "delta_analysis_ternary" / "all_ternary_deltaE_combined.csv",
            output_dir=OUT_ROOT / "ternary",
            required_delta_cols=("chemsys", "source", "delta_E", "metal_1", "metal_2", "ligand"),
        ),
    ]

    for spec in specs:
        print(f"Running {spec.name} analysis...")
        analyze_dataset(spec)
        print(f"Finished {spec.name}. Results: {spec.output_dir}")

    print(f"All results saved to: {OUT_ROOT}")


if __name__ == "__main__":
    main()

