# Maximum-delta_E Point Results

This folder keeps only the single numerically largest `delta_E` row in each category.

Interpretation used here:

- binary category: `(metal, ligand)`
- ternary category: `(metal pair, ligand)`
- `maximum delta_E`: the least-negative or closest-to-zero retained value

## Why this view is useful

These plots answer a best-case question: if each system is represented only by its most favorable available point, how negative can it still remain?

That comparison is especially useful for separating systems that merely have deep excursions from systems that remain persistently negative even in their best-case state.

## Dataset size

- binary selected rows: `37`
- ternary selected rows: `198`

## Files

- `tables/binary_least_negative_points.csv`
- `tables/ternary_least_negative_points.csv`
- `plots/binary_max_delta_points.png`
- `plots/ternary_max_delta_points.png`
- `plots/binary_max_delta_points.html`
- `plots/ternary_max_delta_points.html`

