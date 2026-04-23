# Binary Electronegativity Discussion

This note summarizes the cleaner electronegativity regressions for the binary oxide and sulfide systems.

## Figures used

- `plots/binary_o_clean_multiple_regression.png`
- `plots/binary_s_clean_multiple_regression.png`
- `plots/binary_clean_statistic_comparison_grid.png`
- `tables/layered_regression_summary.csv`

## Main regression summary

| Ligand | Statistic | Correlation (r) | Slope |
|---|---|---:|---:|
| O | min | 0.315 | 0.3384 |
| O | max | -0.258 | -0.0878 |
| O | mean | 0.162 | 0.0463 |
| O | median | 0.126 | 0.0378 |
| S | min | 0.206 | 0.1813 |
| S | max | -0.515 | -0.0990 |
| S | mean | 0.003 | 0.0008 |
| S | median | -0.027 | -0.0077 |

## Interpretation

The clearest trends appear in the extreme statistics rather than in the center of the distribution.

- For oxygen, the most visible trend is the positive slope in `min delta_E`.
- For sulfur, the strongest relation is the negative slope in `max delta_E`.
- Mean and median remain weak for both ligands.

This pattern suggests that electronegativity influences the limiting behavior of the binary `delta_E` distribution more than the average response.

## Practical conclusion

Electronegativity is informative, but it is not the main organizing variable in this project. It should be discussed alongside structural descriptors such as packing fraction rather than used as a standalone explanation.

