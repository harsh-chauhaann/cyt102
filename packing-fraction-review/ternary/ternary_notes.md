# Ternary Packing Efficiency vs delta_E

This folder extends the packing-efficiency analysis to ternary metal-pair systems.

## Dataset summary

- rows in merged table: `1186`
- unique chemsys: `198`
- unique configurations: `1186`
- overall corr(packing_fraction, delta_E): `-0.4795`
- chemsys with usable within-chemsys correlation: `121`
- parse or matching errors: `0`

## Core tables

- `structure_packing_metrics.csv`: one row per configuration with packing descriptors
- `packing_delta_joined.csv`: merged packing and `delta_E` data
- `system_level_overview.csv`: one row per chemsys for inter-compound comparison
- `packing_match_issues.csv`: parse or matching failures

## Reading the results

- use `packing_delta_joined.csv` for configuration-level analysis
- use `system_level_overview.csv` for equal-weight chemsys comparison
- use `packing_vs_delta_corr` and `packing_vs_delta_slope` for within-chemsys trends

## Strongest negative within-chemsys correlations

Many of the most extreme ternary correlations are exactly `-1.0000`, but most of those cases contain only `2` configurations and should therefore be interpreted cautiously.

- `Cd-Mn-S`: `-1.0000` with `2` configurations
- `Hf-Zr-O`: `-1.0000` with `2` configurations
- `Re-Sc-O`: `-1.0000` with `2` configurations
- `Ti-Y-S`: `-1.0000` with `2` configurations
- `Cd-V-O`: `-1.0000` with `2` configurations
- `Ni-Rh-O`: `-1.0000` with `2` configurations
- `Nb-Zr-O`: `-1.0000` with `2` configurations
- `Hf-Mo-S`: `-1.0000` with `2` configurations

## Strongest positive within-chemsys correlations

The same caution applies on the positive side, where several `+1.0000` values also come from two-point systems.

- `Cd-Ti-O`: `1.0000` with `2` configurations
- `Rh-Sc-O`: `1.0000` with `2` configurations
- `Co-Ta-S`: `1.0000` with `2` configurations
- `Cr-Nb-S`: `1.0000` with `2` configurations
- `Nb-Ti-S`: `1.0000` with `2` configurations
- `Fe-Ta-S`: `1.0000` with `2` configurations
- `Mo-Zn-O`: `1.0000` with `2` configurations
- `Co-Mo-O`: `1.0000` with `2` configurations

## Interpretation in one line

Packing fraction shows an even stronger global relationship with ternary `delta_E` than in the binary set, but the ternary space also contains more persistent high-magnitude outliers and more small-sample edge cases.

