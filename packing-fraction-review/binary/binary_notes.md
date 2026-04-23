# Binary Packing Efficiency vs delta_E

This folder links binary `delta_E` behavior to a packing-fraction proxy derived from covalent radii and cell volume.

## Dataset summary

- rows in merged table: `848`
- unique chemsys: `37`
- unique configurations: `848`
- overall corr(packing_fraction, delta_E): `-0.4156`
- chemsys with usable within-chemsys correlation: `32`
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

- `Mn-S`: `-0.8895` with `5` configurations
- `V-S`: `-0.8381` with `13` configurations
- `Sc-O`: `-0.8347` with `8` configurations
- `Fe-S`: `-0.8190` with `12` configurations
- `Y-S`: `-0.8011` with `4` configurations
- `Ta-S`: `-0.7823` with `6` configurations
- `Cr-S`: `-0.7185` with `7` configurations
- `Ti-S`: `-0.7046` with `14` configurations
- `Ta-O`: `-0.6899` with `6` configurations
- `Nb-S`: `-0.6390` with `13` configurations

## Strongest positive within-chemsys correlations

- `Zn-O`: `0.9709` with `5` configurations
- `Hf-O`: `0.8925` with `8` configurations
- `Zn-S`: `0.7103` with `98` configurations
- `Pd-O`: `0.6167` with `3` configurations
- `Pd-S`: `0.3166` with `3` configurations
- `Ni-S`: `0.2586` with `12` configurations
- `Nb-O`: `0.1876` with `17` configurations
- `Sc-S`: `0.1622` with `3` configurations

## Interpretation in one line

Packing fraction has a meaningful global relationship with binary `delta_E`, but it is not a sufficient descriptor on its own because several systems remain strongly negative even at moderate or high packing.

