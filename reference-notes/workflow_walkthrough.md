# Workflow Guide

This project can reproduce the workflow from Materials Project structures to the final binary and ternary summary plots.

## Data source

- Website: <https://materialsproject.org/>
- API documentation: <https://materialsproject.org/api>
- Required credential: `MP_API_KEY`

## Python environment

Install the required Python packages with:

```bash
pip install -r python_packages.txt
```

The main libraries used in this project are:

- `mp-api`
- `pymatgen`
- `pandas`
- `numpy`
- `matplotlib`
- `plotly`
- `adjustText`

## External executable

The workflow also requires Grimme DFT-D3 through the `s-dftd3` command.

Set its location with:

```bash
S_DFTD3_CMD=/path/to/s-dftd3
```

If the executable is already on `PATH`, the default `S_DFTD3_CMD=s-dftd3` is sufficient.

## Configuration

The shared configuration loader reads environment variables directly and also loads a local `.env` file when one is present.

Start from:

```bash
cp sample_environment.env .env
```

Important settings include:

- `MP_API_KEY`
- `S_DFTD3_CMD`
- `MP_TASK_WORKERS`
- `DFTD3_WORKERS`
- `DFTD3_TIMEOUT_SECONDS`
- `MP_API_SLEEP_SECONDS`
- `STRAIN_START`
- `STRAIN_STOP`
- `PIPELINE_METALS`
- `PIPELINE_LIGANDS`
- `ALLOW_IDENTICAL_METAL_PAIRS`

## Full workflow

Run the complete pipeline with:

```bash
python workflow-tools/execute_everything.py
```

Useful variants:

```bash
python workflow-tools/execute_everything.py --binary-only
python workflow-tools/execute_everything.py --ternary-only
python workflow-tools/execute_everything.py --skip-interactive-plots
python workflow-tools/execute_everything.py --dry-run
```

## Step-by-step execution

Binary workflow:

```bash
python binary-structure-sweep/collect_binary_candidates.py
python binary-structure-sweep/measure_binary_dips.py
python binary-structure-sweep/render_binary_rankings.py
```

Ternary workflow:

```bash
python ternary-structure-sweep/collect_ternary_candidates.py
python ternary-structure-sweep/measure_ternary_dips.py
python ternary-structure-sweep/render_ternary_rankings.py
```

Summary analyses:

```bash
python workflow-tools/study_builders/assemble_packing_review.py
python workflow-tools/study_builders/collect_least_negative_snapshots.py
python workflow-tools/study_builders/assemble_aggregate_bars.py
python workflow-tools/study_builders/build_pauling_trend_checks.py
python workflow-tools/study_builders/build_layered_pauling_regressions.py
python workflow-tools/focus_plot_generators/make_binary_focus_panels.py --all
python workflow-tools/focus_plot_generators/make_binary_focus_interactive.py --all
python workflow-tools/focus_plot_generators/make_ternary_focus_interactive.py --all
```

## Output locations

- `binary-structure-sweep/`
- `ternary-structure-sweep/`
- `packing-fraction-review/`
- `least-negative-delta-snapshots/`
- `aggregate-comparison-bars/`
- `pauling-trend-checks/`
- `layered-pauling-regressions/`
- `result_interpretation.md`

