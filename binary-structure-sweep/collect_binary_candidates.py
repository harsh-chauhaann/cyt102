#!/usr/bin/env python3
"""
Binary batch pipeline for fetching Materials Project data and generating delta plots.
"""

import os
import sys
import time
import csv
import re
import shutil
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow-tools"))

from runtime_settings import env_csv_list, env_float, env_int, get_api_key, get_dftd3_command

# ----------------- External libraries -----------------
try:
    from mp_api.client import MPRester
    from pymatgen.io.cif import CifWriter
    from pymatgen.core import Structure
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as e:
    print("Missing dependency:", e)
    print("Install requirements: pip install mp-api pymatgen pandas matplotlib")
    sys.exit(1)

# ======================================================
#                  HARD-CODED CONFIG
# ======================================================

# ---- Metals + Ligands ----
DEFAULT_D_BLOCK_METALS = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",   # 3d
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",   # 4d
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"          # 5d
]
D_BLOCK_METALS = env_csv_list("PIPELINE_METALS", DEFAULT_D_BLOCK_METALS)
LIGANDS = env_csv_list("PIPELINE_LIGANDS", ["O", "S"])

# ---- Directories ----
ROOT_DIR = Path(__file__).resolve().parent
CSV_DIR = ROOT_DIR / "csv"
POSCAR_DIR = ROOT_DIR / "poscars"
STRAIN_DIR = ROOT_DIR / "strains"
AGG_DIR = ROOT_DIR / "aggregated"
PLOT_DIR = ROOT_DIR / "plots"

# ---- DFTD3 ----
S_DFTD3_CMD = get_dftd3_command()
TIMEOUT_RUN = env_int("DFTD3_TIMEOUT_SECONDS", 900)

# ---- Strain grid ----
STRAIN_START = env_int("STRAIN_START", -5)
STRAIN_STOP = env_int("STRAIN_STOP", 29)
STRAIN_INDICES = list(range(STRAIN_START, STRAIN_STOP + 1))

# ---- MP API throttle ----
SLEEP_BETWEEN_API = env_float("MP_API_SLEEP_SECONDS", 0.25)

# ---- Parallelism tuned for 40 cores ----
# Outer: number of concurrent chemsys tasks (metal×ligand)
N_TASK_WORKERS = env_int("MP_TASK_WORKERS", 8)
# Inner: number of concurrent s-dftd3 runs per POSCAR
N_DFTD3_WORKERS = env_int("DFTD3_WORKERS", 4)

# ---- MP fields ----
FIELDS = [
    "material_id", "formula_pretty", "symmetry.crystal_system",
    "symmetry.symbol", "symmetry.number", "nsites",
    "energy_above_hull", "formation_energy_per_atom", "is_stable",
    "volume", "density", "band_gap", "is_gap_direct", "is_metal",
    "ordering", "total_magnetization",
]

# ======================================================
#                     HELPERS
# ======================================================

def ensure_dirs():
    for d in [ROOT_DIR, CSV_DIR, POSCAR_DIR, STRAIN_DIR, AGG_DIR, PLOT_DIR]:
        d.mkdir(parents=True, exist_ok=True)

def sanitize_filename(s: str) -> str:
    if s is None:
        return "unknown"
    s = str(s).strip()
    s = re.sub(r"[^A-Za-z0-9._-]", "_", s)
    return s

# ======================================================
#               MATERIALS PROJECT FETCH
# ======================================================

def fetch_chemsys_write_csv(chemsys: str, api_key: str) -> Path:
    """
    Fetch summary docs and structures from MP for chemsys and write to CSV.
    """
    out_csv = CSV_DIR / f"{chemsys.replace(' ', '')}_table_export_with_structures.csv"

    print(f"[{chemsys}] Fetching from Materials Project...")
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(chemsys=chemsys, fields=FIELDS, all_fields=False)
        docs = list(docs)

    print(f"[{chemsys}] Found {len(docs)} materials")

    rows = []
    for i, d in enumerate(docs, 1):
        mid = d.material_id

        try:
            with MPRester(api_key) as mpr:
                struct = mpr.materials.get_structure_by_material_id(mid)
            cif_text = str(CifWriter(struct))
        except Exception as exc:
            print(f"[{chemsys}] Warning: can't fetch structure {mid}: {exc}")
            cif_text = ""

        row = {
            "Synthesizable": True,
            "Material ID": mid,
            "Formula": d.formula_pretty,
            "Crystal System": getattr(d.symmetry, "crystal_system", None),
            "Space Group Symbol": getattr(d.symmetry, "symbol", None),
            "Space Group Number": getattr(d.symmetry, "number", None),
            "Sites": d.nsites,
            "Energy Above Hull": d.energy_above_hull,
            "Formation Energy": d.formation_energy_per_atom,
            "Predicted Stable": d.is_stable,
            "Volume": d.volume,
            "Density": d.density,
            "Band Gap": d.band_gap,
            "Is Gap Direct": d.is_gap_direct,
            "Is Metal": d.is_metal,
            "Magnetic Ordering": d.ordering,
            "Total Magnetization": d.total_magnetization,
            "Structure": cif_text,
        }
        rows.append(row)

        if i % 10 == 0:
            print(f"[{chemsys}] {i}/{len(docs)} fetched...")

        time.sleep(SLEEP_BETWEEN_API)

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False, quoting=csv.QUOTE_ALL)
    print(f"[{chemsys}] Wrote CSV: {out_csv}")
    return out_csv

# ======================================================
#                 CSV -> POSCAR
# ======================================================

def csv_to_poscars(csv_file: Path, out_poscar_dir: Path) -> int:
    out_poscar_dir.mkdir(parents=True, exist_ok=True)
    written = 0

    with open(csv_file, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SystemExit("No header found in CSV.")

        low_to_field = {
            name.lower().replace(" ", "").replace("_", ""): name
            for name in reader.fieldnames
        }

        def find_field(candidates):
            for cand in candidates:
                key = cand.lower().replace(" ", "").replace("_", "")
                if key in low_to_field:
                    return low_to_field[key]
            for name in reader.fieldnames:
                lname = name.lower()
                for cand in candidates:
                    if cand.lower() in lname:
                        return name
            return None

        formula_key = find_field(["formula"])
        structure_key = find_field(["structure", "cif", "structure_str"])
        material_key = find_field(["materialid", "material id", "material", "mpid"])

        if not formula_key or not structure_key:
            raise SystemExit("Couldn't find required columns in CSV.")

        for idx, row in enumerate(reader):
            try:
                formula = row.get(formula_key, "").strip() or f"unknown_{idx}"
                structure_str = row.get(structure_key, "")
                if not structure_str:
                    continue

                mid_raw = row.get(material_key, "").strip() if material_key else ""
                mid = mid_raw.replace("mp-", "").strip() if mid_raw else str(idx)

                safe_formula = sanitize_filename(formula)
                safe_mid = sanitize_filename(mid)

                cif_name = f"{safe_formula}_{safe_mid}.cif"
                cif_path = out_poscar_dir / cif_name

                with open(cif_path, "w", encoding="utf-8") as cf:
                    cf.write(structure_str)

                try:
                    structure = Structure.from_file(str(cif_path))
                except Exception:
                    cif_path.unlink(missing_ok=True)
                    continue

                poscar_name = f"{safe_formula}_{safe_mid}_POSCAR"
                poscar_path = out_poscar_dir / poscar_name

                with open(poscar_path, "w", encoding="utf-8") as pf:
                    pf.write(structure.to(fmt="poscar"))

                written += 1

            except Exception as exc:
                print(f"[poscar] Failed for row {idx}: {exc}")

    print(f"[poscar] Total POSCARs created: {written}")
    return written

# ======================================================
#            s-dftd3 OUTPUT PARSING
# ======================================================

def parse_output_for_props(
    text: str,
    cation_symbol: Optional[str] = None,
    anion_symbol: Optional[str] = None
) -> Tuple[str, str, str, str, str]:
    """
    Extract CN/C6 for cation and anion and the dispersion energy.
    Returns: CN_c, CN_a, C6_c, C6_a, Edis
    """
    CN_c = CN_a = C6_c = C6_a = "N/A"
    Edis = "N/A"
    lines = text.splitlines()

    # Dispersion energy
    for ln in lines:
        m = re.search(r"Dispersion\s+energy[:\s]+([-\d\.Ee+]+)", ln, re.I)
        if m:
            Edis = m.group(1)
            break

    # Property table parse (your original logic)
    elem_map = {}
    for ln in lines:
        toks = ln.split()
        if len(toks) >= 5:
            sym = toks[2]
            cn = toks[3]
            c6 = toks[4]
            if sym not in elem_map:
                elem_map[sym] = (cn, c6)

    # Cation
    if cation_symbol and cation_symbol in elem_map:
        CN_c, C6_c = elem_map[cation_symbol]
    elif elem_map:
        first = next(iter(elem_map))
        CN_c, C6_c = elem_map[first]

    # Anion
    if anion_symbol and anion_symbol in elem_map:
        CN_a, C6_a = elem_map[anion_symbol]
    elif len(elem_map) >= 2:
        second = list(elem_map.keys())[1]
        CN_a, C6_a = elem_map[second]
    elif elem_map:
        CN_a, C6_a = elem_map[next(iter(elem_map))]

    return CN_c, CN_a, C6_c, C6_a, Edis

# ======================================================
#             STRAIN + s-dftd3 EXECUTION
# ======================================================

def run_single_strain(args):
    """
    Run s-dftd3 for a single strain index.
    This function is picklable and safe for ProcessPoolExecutor.
    """
    i, base_vectors, orig_lines, target_dir, results_dir, s_dftd3_cmd = args

    scale = 0.05 * i
    factor = 1.0 + scale
    new_vectors = [[v * factor for v in vec] for vec in base_vectors]

    poscar_name = f"POSCAR{i}"
    poscar_path = target_dir / poscar_name

    new_lines = list(orig_lines)
    while len(new_lines) < 5:
        new_lines.append("")

    def fmt_vec(v):
        return "{:18.10f} {:18.10f} {:18.10f}".format(v[0], v[1], v[2])

    new_lines[2] = fmt_vec(new_vectors[0])
    new_lines[3] = fmt_vec(new_vectors[1])
    new_lines[4] = fmt_vec(new_vectors[2])

    with open(poscar_path, "w") as wf:
        wf.write("\n".join(new_lines) + "\n")

    out_path = results_dir / f"out{i}.txt"
    cmd = [
        s_dftd3_cmd,
        "-i", "vasp",
        str(poscar_path),
        "--zero", "pbe",
        "--property",
        "--verbose",
    ]

    try:
        if shutil.which(s_dftd3_cmd) is None and not os.path.exists(s_dftd3_cmd):
            with open(out_path, "w") as outf:
                outf.write(f"=== RUN SKIPPED: {s_dftd3_cmd} not found in PATH ===\n")
        else:
            with open(out_path, "w") as outf:
                subprocess.run(cmd, stdout=outf, stderr=subprocess.STDOUT, timeout=TIMEOUT_RUN)
    except Exception as e:
        with open(out_path, "a") as outf:
            outf.write(f"\n\n=== RUN ERROR ===\n{e}\n")

    return i, poscar_path, out_path

def process_single_poscar(poscar_file: Path, chemsys: str, n_workers: int) -> List[dict]:
    """
    For one POSCAR, generate strained POSCARs and run s-dftd3 for each strain in parallel.
    """
    base_name = poscar_file.stem
    target_dir = STRAIN_DIR / chemsys.replace(" ", "") / base_name
    results_dir = target_dir / "results"
    target_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    try:
        with open(poscar_file, "r", errors="ignore") as pf:
            orig_lines = [ln.rstrip("\n") for ln in pf.readlines()]
    except Exception as e:
        print(f"[{chemsys}] Can't read {poscar_file}: {e}")
        return []

    if len(orig_lines) < 5:
        return []

    def parse_vector_line(line):
        parts = line.split()
        vals = []
        for p in parts[:3]:
            try:
                vals.append(float(p))
            except:
                vals.append(0.0)
        while len(vals) < 3:
            vals.append(0.0)
        return vals

    base_vectors = [
        parse_vector_line(orig_lines[2]),
        parse_vector_line(orig_lines[3]),
        parse_vector_line(orig_lines[4]),
    ]

    # element list is line 6 in VASP POSCAR format
    elements = []
    if len(orig_lines) >= 6:
        parts = orig_lines[5].strip().split()
        if parts:
            elements = parts

    cation_symbol = elements[0] if elements else None
    anion_symbol = elements[-1] if elements else None

    strain_args = [
        (i, base_vectors, orig_lines, target_dir, results_dir, S_DFTD3_CMD)
        for i in STRAIN_INDICES
    ]

    print(f"[{chemsys}] Processing {poscar_file.name} with {n_workers} workers...")
    results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(run_single_strain, arg) for arg in strain_args]
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                print(f"[{chemsys}] Strain calculation failed: {e}")

    rows = []
    for i, poscar_path, out_path in results:
        # cell_length from POSCAR line 3 (index 2), first number
        try:
            with open(poscar_path, "r") as pf:
                p_lines = pf.readlines()
            cell_length = p_lines[2].split()[0] if len(p_lines) >= 3 else "N/A"
        except:
            cell_length = "N/A"

        try:
            out_text = open(out_path, "r", errors="ignore").read()
        except:
            out_text = ""

        CN_c, CN_a, C6_c, C6_a, Edis = parse_output_for_props(
            out_text, cation_symbol, anion_symbol
        )

        rows.append({
            "source": base_name,
            "cell_length": cell_length,
            "CN_C": CN_c,
            "CN_A": CN_a,
            "C6_C": C6_c,
            "C6_A": C6_a,
            "Edis": Edis,
        })

    return rows

def run_strain_and_collect(poscar_dir: Path, chemsys: str) -> Path:
    """
    Process all POSCARs sequentially, but each POSCAR runs strains in parallel.
    """
    agg_csv = AGG_DIR / f"{chemsys.replace(' ', '')}_output.csv"
    fieldnames = ["source", "cell_length", "CN_C", "CN_A", "C6_C", "C6_A", "Edis"]

    poscar_files = [f for f in sorted(poscar_dir.iterdir()) if f.is_file()]
    all_rows = []

    for poscar_file in poscar_files:
        rows = process_single_poscar(poscar_file, chemsys, N_DFTD3_WORKERS)
        all_rows.extend(rows)

    with open(agg_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"[{chemsys}] Aggregated CSV written: {agg_csv}")
    return agg_csv

# ======================================================
#                     PLOTTING
# ======================================================

def plot_aggregated_csv(agg_csv: Path, chemsys: str, y_min: float = -0.6, y_max: float = 0.0):
    df = pd.read_csv(agg_csv)
    df["cell_length"] = pd.to_numeric(df["cell_length"], errors="coerce")
    df["Edis"] = pd.to_numeric(df["Edis"], errors="coerce")
    df = df.dropna(subset=["cell_length", "Edis", "source"])

    if df.empty:
        print(f"[{chemsys}] No valid data to plot")
        return

    fig, ax = plt.subplots(figsize=(9, 6))
    for source, grp in df.groupby("source"):
        ax.plot(grp["cell_length"], grp["Edis"], marker="o", linestyle="-", label=source)

    ax.set_xlabel("Cell Length")
    ax.set_ylabel("Dispersion Energy (Edis)")
    ax.set_ylim(y_min, y_max)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
    plt.title(f"{chemsys}: Edis vs cell length")
    plt.tight_layout()

    out_png = PLOT_DIR / f"{chemsys.replace(' ', '')}_cell_length_vs_Edis.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[{chemsys}] Saved plot: {out_png}")

# ======================================================
#               PER CHEMSYS PIPELINE
# ======================================================

def process_chemsys(metal: str, ligand: str, api_key: str):
    chemsys = f"{metal}-{ligand}"
    chemsys_tag = chemsys.replace(" ", "")

    print(f"\n{'='*60}")
    print(f"[{chemsys}] START")
    print(f"{'='*60}")

    # Checkpoint: if aggregated results exist, skip
    agg_csv = AGG_DIR / f"{chemsys_tag}_output.csv"
    if agg_csv.exists():
        print(f"[{chemsys}] Found existing aggregated CSV -> SKIP: {agg_csv}")
        return chemsys, "SKIPPED (agg exists)"

    try:
        poscar_out_dir = POSCAR_DIR / chemsys_tag

        # 1) Fetch MP data + structures -> CSV
        csv_path = fetch_chemsys_write_csv(chemsys, api_key)

        # 2) CSV -> POSCAR
        poscar_out_dir.mkdir(parents=True, exist_ok=True)
        csv_to_poscars(csv_path, poscar_out_dir)

        # 3) Strains + s-dftd3 -> aggregated CSV
        agg_csv = run_strain_and_collect(poscar_out_dir, chemsys)

        # 4) Plot
        plot_aggregated_csv(agg_csv, chemsys)

        print(f"[{chemsys}] DONE ✓")
        return chemsys, "SUCCESS"

    except Exception as e:
        print(f"[{chemsys}] ERROR: {e}")
        import traceback
        traceback.print_exc()
        return chemsys, f"FAILED: {e}"

# ======================================================
#                      MAIN
# ======================================================

def main():
    ensure_dirs()

    api_key = get_api_key()
    if not api_key:
        print("Materials Project API key not found.")
        print("Create one at https://materialsproject.org/api and set MP_API_KEY")
        print("or place it in ~/.mp_api_key (or the file pointed to by MP_API_KEY_FILE).")
        sys.exit(1)

    # Print config summary
    import multiprocessing
    print("\n==============================")
    print("HARD-CODED PIPELINE CONFIG")
    print("==============================")
    print(f"Metals (d-block): {len(D_BLOCK_METALS)}")
    print(f"Ligands: {LIGANDS}")
    print(f"Total chemsys tasks: {len(D_BLOCK_METALS) * len(LIGANDS)}")
    print(f"Output directory: {ROOT_DIR}")
    print(f"s-dftd3 cmd: {S_DFTD3_CMD}")
    print(f"Strain indices: {STRAIN_START} .. {STRAIN_STOP}")
    print(f"Available CPUs: {multiprocessing.cpu_count()}")
    print(f"Outer task workers: {N_TASK_WORKERS}")
    print(f"Inner dftd3 workers per POSCAR: {N_DFTD3_WORKERS}")
    print(f"Total potential parallelism: {N_TASK_WORKERS * N_DFTD3_WORKERS}")
    print(f"MP API sleep: {SLEEP_BETWEEN_API} sec")
    print("==============================\n")

    # Build task list (metal×ligand)
    tasks = [(m, l) for m in D_BLOCK_METALS for l in LIGANDS]

    # Run tasks in parallel
    results = []
    with ProcessPoolExecutor(max_workers=N_TASK_WORKERS) as executor:
        futures = {executor.submit(process_chemsys, m, l, api_key): (m, l) for (m, l) in tasks}

        for future in as_completed(futures):
            m, l = futures[future]
            chemsys = f"{m}-{l}"
            try:
                results.append(future.result())
            except Exception as e:
                print(f"[{chemsys}] FAILED with exception: {e}")
                results.append((chemsys, f"FAILED: {e}"))

    # Summary
    print(f"\n{'='*60}")
    print("PIPELINE COMPLETE - Summary:")
    print(f"{'='*60}")
    for chemsys, status in sorted(results, key=lambda x: x[0]):
        print(f"  {chemsys}: {status}")

    print(f"\nResults directory: {ROOT_DIR}\n")

if __name__ == "__main__":
    main()

