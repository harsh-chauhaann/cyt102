#!/usr/bin/env python3
"""
measure_binary_dips.py  —  Dispersion-energy dip analysis
-------------------------------------------------
Reads aggregated *_output.csv files produced by collect_binary_candidates.py,
detects the bond-length dip in each Edis vs cell_length curve,
estimates a linear baseline, and records ΔE = E_min − E_expected.

Dip detection uses TWO strategies:
  1. Local-max-then-drop  — robust against rising/falling baselines
  2. Global gradient threshold — fallback for simple bowl-shaped curves
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import traceback
from collections import defaultdict

# ============================================================
#  PATHS
# ============================================================
ROOT     = Path(__file__).resolve().parent
AGG_DIR  = ROOT / "aggregated"
OUT_DIR  = ROOT / "delta_analysis"
PLOT_DIR = OUT_DIR / "plots"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
#  TUNEABLE PARAMETERS
# ============================================================
MIN_POINTS          = 10     # minimum data points per source curve
SLOPE_SIGMA         = 2.0    # std-deviations below mean gradient (fallback dip detection)
RIGHT_FLAT_PCTILE   = 25     # percentile of |dy| on right side treated as "flat"
RIGHT_MIN_OFFSET    = 2      # minimum points past dip before looking for right anchor
LEFT_FALLBACK_FRAC  = 0.33   # fraction of dip_idx used as fallback left offset
DELTA_E_MAX         = 0.0    # only keep dips with delta_E < this value
MIN_DIP_DROP        = 0.05   # minimum y drop from local peak to dip to count as real
LOCAL_MAX_WINDOW    = 15     # how many points after a peak to search for dip minimum

# ============================================================
#  SKIP-REASON COUNTERS
# ============================================================
skip_counts: dict = defaultdict(int)


def record_skip(reason: str, chemsys: str, source: str, extra: str = ""):
    skip_counts[reason] += 1
    msg = f"  [SKIP:{reason}] {chemsys} | {source}"
    if extra:
        msg += f"  →  {extra}"
    print(msg)


# ============================================================
#  STEP 1 – CLEAN
# ============================================================
def clean_data(grp: pd.DataFrame) -> pd.DataFrame:
    grp = grp.copy()
    grp["cell_length"] = pd.to_numeric(grp["cell_length"], errors="coerce")
    grp["Edis"]        = pd.to_numeric(grp["Edis"],        errors="coerce")
    grp = grp.replace([np.inf, -np.inf], np.nan)
    grp = grp.dropna(subset=["cell_length", "Edis"])
    grp = grp.sort_values("cell_length").reset_index(drop=True)
    return grp


# ============================================================
#  STEP 2 – DIP DETECTION  (two strategies)
# ============================================================
def detect_dip(x: np.ndarray, y: np.ndarray):
    """
    Returns (dip_idx, debug_dict) or (None, debug_dict).

    Strategy 1 — Local max then drop:
        Find every local maximum in y, then search the next
        LOCAL_MAX_WINDOW points for the deepest drop.
        Works even when the whole curve has a strong rising or
        falling background (e.g. Zr64S127 style curves).

    Strategy 2 — Global gradient threshold (fallback):
        Classic: find points where dy << mean(dy) - SLOPE_SIGMA*std(dy).
        Works for simple bowl-shaped curves with flat backgrounds.
    """
    debug = {}
    n     = len(y)
    dy    = np.gradient(y, x)

    # ── Strategy 1: local max → drop ─────────────────────────
    local_maxima = []
    for i in range(1, n - 1):
        if y[i] > y[i - 1] and y[i] > y[i + 1]:
            local_maxima.append(i)

    debug["n_local_maxima"] = len(local_maxima)

    dip_idx   = None
    best_drop = 0.0
    peak_used = None

    for peak_i in local_maxima:
        search_end = min(peak_i + LOCAL_MAX_WINDOW, n)
        for j in range(peak_i + 1, search_end):
            drop = y[peak_i] - y[j]
            if drop > best_drop and drop > MIN_DIP_DROP:
                best_drop = drop
                dip_idx   = j
                peak_used = peak_i

    if dip_idx is not None:
        debug["method"]         = "local_max_drop"
        debug["peak_idx"]       = int(peak_used)
        debug["peak_x"]         = round(float(x[peak_used]), 6)
        debug["peak_y"]         = round(float(y[peak_used]), 6)
        debug["drop_magnitude"] = round(best_drop, 6)
        debug["dip_idx"]        = int(dip_idx)
        debug["dip_x"]          = round(float(x[dip_idx]), 6)
        debug["dip_y"]          = round(float(y[dip_idx]), 6)
        debug["dip_gradient"]   = round(float(dy[dip_idx]), 6)
        return dip_idx, debug

    # ── Strategy 2: global gradient threshold (fallback) ─────
    mean_dy   = float(np.mean(dy))
    std_dy    = float(np.std(dy))
    threshold = mean_dy - SLOPE_SIGMA * std_dy

    debug["method"]       = "global_gradient_fallback"
    debug["dy_mean"]      = round(mean_dy, 6)
    debug["dy_std"]       = round(std_dy, 6)
    debug["threshold"]    = round(threshold, 6)

    candidates = np.where(dy < threshold)[0]
    debug["n_candidates"] = int(len(candidates))

    if len(candidates) == 0:
        debug["reason"] = "no local maxima with sufficient drop AND no gradient candidates"
        return None, debug

    dip_idx = int(candidates[np.argmin(y[candidates])])
    debug["dip_idx"]      = dip_idx
    debug["dip_x"]        = round(float(x[dip_idx]), 6)
    debug["dip_y"]        = round(float(y[dip_idx]), 6)
    debug["dip_gradient"] = round(float(dy[dip_idx]), 6)

    return dip_idx, debug


# ============================================================
#  STEP 3 – BASELINE ESTIMATION
# ============================================================
def estimate_baseline(x: np.ndarray, y: np.ndarray, dip_idx: int, peak_idx: int = None):
    """
    Returns (E_expected, y_baseline, debug_dict) or (None, None, debug_dict).

    Left anchor:
        If a peak_idx is known (from Strategy 1), use it directly —
        it is exactly the shoulder before the dip.
        Otherwise walk left from dip_idx to find the gradient shoulder,
        with a fractional fallback.

    Right anchor:
        Flattest quartile of |dy| to the right of the dip.
        Fallback: last 3 points.
    """
    debug = {}
    try:
        n  = len(x)
        dy = np.gradient(y, x)

        # ── LEFT ANCHOR ──────────────────────────────────────
        if peak_idx is not None:
            left_idx = peak_idx
            debug["left_method"] = "peak from local_max_drop"
        else:
            left_idx = None
            for i in range(dip_idx - 1, 0, -1):
                if abs(dy[i]) < abs(dy[i - 1]):
                    left_idx = i
                else:
                    if left_idx is not None:
                        break

            if left_idx is None or left_idx < 1:
                fallback_offset = max(2, int(dip_idx * LEFT_FALLBACK_FRAC))
                left_idx = max(0, dip_idx - fallback_offset)
                debug["left_method"] = f"fallback (offset={fallback_offset})"
            else:
                debug["left_method"] = "gradient shoulder"

        debug["left_idx"] = int(left_idx)
        debug["left_x"]   = round(float(x[left_idx]), 6)
        debug["left_y"]   = round(float(y[left_idx]), 6)

        # ── RIGHT ANCHOR ─────────────────────────────────────
        right_start      = dip_idx + RIGHT_MIN_OFFSET
        right_candidates = []

        if right_start < n:
            dy_right_abs   = np.abs(dy[right_start:])
            flat_threshold = float(np.percentile(dy_right_abs, RIGHT_FLAT_PCTILE))
            flat_threshold = max(flat_threshold, 1e-8)
            debug["right_flat_threshold"] = round(flat_threshold, 8)

            for i in range(right_start, n):
                if abs(dy[i]) <= flat_threshold:
                    right_candidates.append(i)

        debug["n_right_candidates"] = len(right_candidates)

        if not right_candidates:
            right_idx = max(right_start, n - 3)
            debug["right_method"] = "fallback (last 3 pts)"
        else:
            right_idx = right_candidates[-1]
            debug["right_method"] = "flattest quartile"

        window  = min(3, n - right_idx)
        x_right = float(np.mean(x[right_idx: right_idx + window]))
        y_right = float(np.mean(y[right_idx: right_idx + window]))

        debug["right_idx"] = int(right_idx)
        debug["right_x"]   = round(x_right, 6)
        debug["right_y"]   = round(y_right, 6)

        # ── LINEAR BASELINE ───────────────────────────────────
        dx = x_right - float(x[left_idx])
        if abs(dx) < 1e-10:
            debug["reason"] = "left and right anchors at same x"
            return None, None, debug

        slope     = (y_right - float(y[left_idx])) / dx
        intercept = float(y[left_idx]) - slope * float(x[left_idx])

        E_expected = slope * float(x[dip_idx]) + intercept
        y_baseline = slope * x + intercept

        debug["slope"]      = round(slope, 8)
        debug["intercept"]  = round(intercept, 8)
        debug["E_expected"] = round(E_expected, 8)

        return E_expected, y_baseline, debug

    except Exception as e:
        debug["reason"] = f"exception: {e}"
        return None, None, debug


# ============================================================
#  STEP 4 – PROCESS ONE FILE
# ============================================================
def process_file(csv_path: Path):
    chemsys = csv_path.stem.replace("_output", "")
    print(f"\n{'='*64}")
    print(f"  PROCESSING  {chemsys}   ({csv_path.name})")
    print(f"{'='*64}")

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[ERROR] Cannot read file: {e}")
        skip_counts["file_read_error"] += 1
        return

    if df.empty:
        print("[WARNING] Empty CSV — skipping")
        skip_counts["empty_csv"] += 1
        return

    sources    = df["source"].nunique()
    print(f"  Sources in file: {sources}")

    results    = []
    file_skips = defaultdict(int)

    for source, grp in df.groupby("source"):
        print(f"\n  ── {source} ──")

        # ── clean ─────────────────────────────────────────────
        grp = clean_data(grp)
        n   = len(grp)
        print(f"     Points after cleaning : {n}")

        if n < MIN_POINTS:
            record_skip("too_few_points", chemsys, source,
                        f"need {MIN_POINTS}, have {n}")
            file_skips["too_few_points"] += 1
            continue

        x = grp["cell_length"].values
        y = grp["Edis"].values

        print(f"     x range : {x.min():.4f} → {x.max():.4f}")
        print(f"     y range : {y.min():.6f} → {y.max():.6f}")

        # ── dip detection ─────────────────────────────────────
        dip_idx, dip_dbg = detect_dip(x, y)
        print(f"     Dip detection  : {dip_dbg}")

        if dip_idx is None:
            record_skip("no_dip", chemsys, source,
                        dip_dbg.get("reason", ""))
            file_skips["no_dip"] += 1
            continue

        # pass peak_idx to baseline if Strategy 1 was used
        peak_idx = dip_dbg.get("peak_idx", None)

        # ── baseline ──────────────────────────────────────────
        E_expected, y_baseline, bl_dbg = estimate_baseline(x, y, dip_idx, peak_idx)
        print(f"     Baseline debug : {bl_dbg}")

        if E_expected is None:
            record_skip("baseline_failed", chemsys, source,
                        bl_dbg.get("reason", ""))
            file_skips["baseline_failed"] += 1
            continue

        # ── delta E ───────────────────────────────────────────
        E_min   = float(y[dip_idx])
        delta_E = E_min - E_expected
        print(f"     E_min={E_min:.6f}  E_expected={E_expected:.6f}"
              f"  ΔE={delta_E:.6f}")

        if delta_E >= DELTA_E_MAX:
            record_skip("delta_E_positive", chemsys, source,
                        f"delta_E={delta_E:.6f} >= {DELTA_E_MAX}")
            file_skips["delta_E_positive"] += 1
            continue

        print(f"     ✓  ACCEPTED  ΔE = {delta_E:.6f}")

        results.append({
            "chemsys":         chemsys,
            "source":          source,
            "n_points":        n,
            "dip_method":      dip_dbg.get("method", ""),
            "dip_cell_length": float(x[dip_idx]),
            "E_min":           E_min,
            "E_expected":      E_expected,
            "delta_E":         delta_E,
            "left_x":          bl_dbg.get("left_x"),
            "right_x":         bl_dbg.get("right_x"),
            "baseline_slope":  bl_dbg.get("slope"),
        })

        # ── plot ──────────────────────────────────────────────
        try:
            fig, ax = plt.subplots(figsize=(8, 5))

            ax.plot(x, y, "o-", color="steelblue", lw=1.5,
                    ms=4, label="Edis")
            ax.plot(x, y_baseline, "--", color="darkorange", lw=1.5,
                    label="Baseline")

            ax.axvline(x[bl_dbg["left_idx"]], color="green",
                       ls=":", lw=1.2,
                       label=f"left  x={bl_dbg['left_x']:.3f}  [{bl_dbg['left_method']}]")

            ax.axvline(bl_dbg["right_x"], color="purple",
                       ls=":", lw=1.2,
                       label=f"right x={bl_dbg['right_x']:.3f}  [{bl_dbg['right_method']}]")

            if peak_idx is not None:
                ax.scatter(x[peak_idx], y[peak_idx], color="limegreen",
                           s=80, marker="^", zorder=6,
                           label=f"Peak  x={x[peak_idx]:.3f}")

            ax.scatter(x[dip_idx], E_min, color="red", s=90,
                       zorder=7, label=f"Dip  ΔE={delta_E:.4f}")

            ax.scatter(x[dip_idx], E_expected, color="orange",
                       s=70, marker="^", zorder=7, label="E_expected")

            ax.set_xlabel("Cell Length")
            ax.set_ylabel("Edis")
            ax.set_title(f"{chemsys} | {source}\n"
                         f"method={dip_dbg.get('method','')}  "
                         f"ΔE={delta_E:.4f}")
            ax.legend(fontsize=7, loc="best")
            fig.tight_layout()

            out_plot = PLOT_DIR / f"{chemsys}_{source}.png"
            fig.savefig(out_plot, dpi=200, bbox_inches="tight")
            plt.close(fig)
            print(f"     Plot saved : {out_plot.name}")

        except Exception as e:
            print(f"     [WARNING] Plot failed: {e}")
            traceback.print_exc()

    # ── per-file summary ──────────────────────────────────────
    print(f"\n  ── File summary for {chemsys} ──")
    print(f"     Total sources : {sources}")
    print(f"     Accepted      : {len(results)}")
    for reason, cnt in sorted(file_skips.items()):
        print(f"     Skipped [{reason}] : {cnt}")

    # ── save CSV ──────────────────────────────────────────────
    if results:
        out_csv = OUT_DIR / f"{chemsys}_deltaE.csv"
        pd.DataFrame(results).to_csv(out_csv, index=False)
        print(f"  → Saved: {out_csv}")
    else:
        print("  → No valid dips found — no CSV written")


# ============================================================
#  MAIN
# ============================================================
def main():
    files = sorted(AGG_DIR.glob("*_output.csv"))

    if not files:
        print(f"[ERROR] No aggregated CSV files found in {AGG_DIR}")
        return

    print(f"Found {len(files)} aggregated file(s) in {AGG_DIR}\n")

    all_result_rows = []

    for f in files:
        process_file(f)

        chemsys   = f.stem.replace("_output", "")
        delta_csv = OUT_DIR / f"{chemsys}_deltaE.csv"
        if delta_csv.exists():
            try:
                all_result_rows.append(pd.read_csv(delta_csv))
            except Exception:
                pass

    # ── global summary ────────────────────────────────────────
    print(f"\n{'='*64}")
    print("  GLOBAL SKIP SUMMARY")
    print(f"{'='*64}")
    total_skipped = sum(skip_counts.values())
    for reason, cnt in sorted(skip_counts.items(), key=lambda x: -x[1]):
        print(f"  {reason:<30} : {cnt}")
    print(f"  {'TOTAL skipped':<30} : {total_skipped}")

    if all_result_rows:
        combined = pd.concat(all_result_rows, ignore_index=True)
        print(f"\n  Total accepted dips        : {len(combined)}")
        print(f"  Unique chemsys             : {combined['chemsys'].nunique()}")
        print(f"\n  Dip method breakdown:")
        print(combined["dip_method"].value_counts().to_string())
        print(f"\n  ΔE statistics:")
        print(combined["delta_E"].describe().to_string())

        combined_out = OUT_DIR / "all_deltaE_combined.csv"
        combined.to_csv(combined_out, index=False)
        print(f"\n  Combined CSV → {combined_out}")
    else:
        print("\n  No accepted dips across all files.")

    print(f"\n{'='*64}")
    print("  ALL DONE")
    print(f"{'='*64}\n")


if __name__ == "__main__":
    main()

