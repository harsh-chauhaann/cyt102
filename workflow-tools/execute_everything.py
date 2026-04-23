#!/usr/bin/env python3
"""
Run the full fetch-to-plot pipeline from Materials Project through summary plots.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from runtime_settings import REPO_ROOT


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the full binary + ternary delta-analysis pipeline.",
    )
    parser.add_argument("--binary-only", action="store_true", help="Run only the binary pipeline.")
    parser.add_argument("--ternary-only", action="store_true", help="Run only the ternary pipeline.")
    parser.add_argument("--skip-extremes", action="store_true", help="Skip top-vs-bottom figure scripts.")
    parser.add_argument("--skip-focused-plots", action="store_true", help="Skip per-metal static packing-vs-delta plots.")
    parser.add_argument("--skip-interactive-plots", action="store_true", help="Skip Plotly HTML plot generation.")
    parser.add_argument("--dry-run", action="store_true", help="Print the steps without running them.")
    return parser.parse_args()


def script_step(script: str, *extra_args: str) -> list[str]:
    return [sys.executable, str(REPO_ROOT / script), *extra_args]


def run_step(command: list[str], dry_run: bool) -> None:
    print(f"\n>>> {' '.join(command)}")
    if dry_run:
        return
    subprocess.run(command, check=True, cwd=REPO_ROOT)


def build_steps(args: argparse.Namespace) -> list[list[str]]:
    include_binary = not args.ternary_only
    include_ternary = not args.binary_only

    steps: list[list[str]] = []

    if include_binary:
        steps.extend([
            script_step("binary-structure-sweep/collect_binary_candidates.py"),
            script_step("binary-structure-sweep/measure_binary_dips.py"),
        ])
        if not args.skip_extremes:
            steps.append(script_step("binary-structure-sweep/render_binary_rankings.py"))

    if include_ternary:
        steps.extend([
            script_step("ternary-structure-sweep/collect_ternary_candidates.py"),
            script_step("ternary-structure-sweep/measure_ternary_dips.py"),
        ])
        if not args.skip_extremes:
            steps.append(script_step("ternary-structure-sweep/render_ternary_rankings.py"))

    steps.extend([
        script_step("workflow-tools/study_builders/assemble_packing_review.py"),
        script_step("workflow-tools/study_builders/collect_least_negative_snapshots.py"),
        script_step("workflow-tools/study_builders/assemble_aggregate_bars.py"),
        script_step("workflow-tools/study_builders/build_pauling_trend_checks.py"),
    ])

    if not args.skip_focused_plots:
        steps.append(script_step("workflow-tools/focus_plot_generators/make_binary_focus_panels.py", "--all"))

    if not args.skip_interactive_plots:
        steps.append(script_step("workflow-tools/focus_plot_generators/make_binary_focus_interactive.py", "--all"))
        steps.append(script_step("workflow-tools/focus_plot_generators/make_ternary_focus_interactive.py", "--all"))

    return steps


def main() -> None:
    args = parse_args()
    if args.binary_only and args.ternary_only:
        raise SystemExit("Choose only one of --binary-only or --ternary-only.")

    steps = build_steps(args)
    print(f"Pipeline root: {REPO_ROOT}")
    print(f"Total steps: {len(steps)}")

    for command in steps:
        run_step(command, dry_run=args.dry_run)

    if args.dry_run:
        print("\nDry run complete.")
    else:
        print("\nPipeline complete.")


if __name__ == "__main__":
    main()

