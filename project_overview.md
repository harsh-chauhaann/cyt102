# D-Metal Oxide and Sulfide Delta Analysis

This folder is a clean working copy of the delta-analysis project for binary and ternary d-metal oxide and sulfide systems. It combines the reproducible workflow, the processed outputs, and the interpretation notes in one place.

The project has two main goals:

- reproduce the strain-based `delta_E` workflow from Materials Project structures
- summarize the chemical and structural trends visible in the final plots

## Project structure

- `binary-structure-sweep/`: binary fetch, strain, dip-detection, and binary summary outputs
- `ternary-structure-sweep/`: ternary fetch, strain, dip-detection, and ternary summary outputs
- `packing-fraction-review/`: packing-fraction analysis for binary and ternary systems
- `least-negative-delta-snapshots/`: best-case `delta_E` point retained for each binary metal or ternary metal pair
- `aggregate-comparison-bars/`: grouped `min`, `max`, `mean`, and `median` comparisons
- `pauling-trend-checks/`: single-regression electronegativity plots
- `layered-pauling-regressions/`: cleaner multi-regression follow-up figures
- `result_interpretation.md`: main interpretation note across the major figure families
- `reference-notes/workflow_walkthrough.md`: workflow and reproducibility guide

## Data source and requirements

The structure data come from the [Materials Project](https://materialsproject.org/). Running the full workflow requires:

- a valid `MP_API_KEY`
- the `s-dftd3` executable
- the Python packages listed in `python_packages.txt`

Start from:

```bash
pip install -r python_packages.txt
```

Then create a local environment file from:

```bash
cp sample_environment.env .env
```

## Running the workflow

The complete run can be launched with:

```bash
python workflow-tools/execute_everything.py
```

Helpful variants:

```bash
python workflow-tools/execute_everything.py --binary-only
python workflow-tools/execute_everything.py --ternary-only
python workflow-tools/execute_everything.py --skip-interactive-plots
python workflow-tools/execute_everything.py --dry-run
```

## Recommended reading order

If the goal is interpretation rather than rerunning calculations, the most useful order is:

1. `result_interpretation.md`
2. `packing-fraction-review/`
3. `least-negative-delta-snapshots/`
4. `aggregate-comparison-bars/`
5. `pauling-trend-checks/`
6. `layered-pauling-regressions/`

## Summary

The strongest recurring picture in this dataset is that packing efficiency shows a clearer global relationship with `delta_E` than electronegativity does, but neither descriptor fully explains the most negative outliers. Binary systems often approach near-zero best-case values, whereas ternary systems retain more severe and more persistent negative behavior.

