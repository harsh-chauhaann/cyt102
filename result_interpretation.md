# Delta Analysis of Binary and Ternary d-Metal Systems

This note brings together the main interpretations from the major figure sets in this project. The emphasis is on what the plots show about `delta_E`, rather than on the mechanics of how each file was generated.

The discussion below uses the following folders:

- `packing-fraction-review/`
- `least-negative-delta-snapshots/`
- `aggregate-comparison-bars/`
- `pauling-trend-checks/`
- `layered-pauling-regressions/`

## 1. Packing efficiency and delta_E

### Binary figures

![Binary configuration-level packing vs delta_E](packing-fraction-review/binary/plots/configuration_level_packing_vs_delta.png)

![Binary inter-compound packing comparison](packing-fraction-review/binary/plots/chemsys_level_median_packing_vs_best_delta.png)

![Binary within-chemsys correlation histogram](packing-fraction-review/binary/plots/within_chemsys_correlation_histogram.png)

### Ternary figures

![Ternary configuration-level packing vs delta_E](packing-fraction-review/ternary/plots/configuration_level_packing_vs_delta.png)

![Ternary inter-compound packing comparison](packing-fraction-review/ternary/plots/chemsys_level_median_packing_vs_best_delta.png)

![Ternary within-chemsys correlation histogram](packing-fraction-review/ternary/plots/within_chemsys_correlation_histogram.png)

### Interpretation

The binary configuration-level scatter shows a broad negative tendency: many data points move below zero as packing fraction increases, and some of the deepest negative outliers occur in the medium-to-high packing range. This means packing fraction matters, but the plot does not support the simplified idea that denser packing automatically suppresses the `delta_E` dip.

The binary inter-compound comparison makes that point even more clearly. `Ti-O` is the strongest outlier, reaching roughly `-1.61`, yet it is not a low-packing system. Other labeled systems such as `Nb-S`, `Ti-S`, `Zr-S`, `V-O`, `Y-O`, and `Cu-O` also remain strongly negative across a wide span of median packing fractions. The figure therefore suggests that chemistry and local motif are still crucial even when packing is relatively high.

The binary within-chemsys histogram is weighted toward negative correlations, which supports the general pattern that higher packing often corresponds to lower `delta_E` within a given binary system. At the same time, the positive side is not empty. A few systems clearly behave differently, so packing fraction should be treated as an informative descriptor rather than a universal control parameter.

The ternary configuration-level plot follows the same overall direction, but more strongly. Oxygen ternaries extend further into the deep negative region, while sulfur ternaries are somewhat more compressed but still clearly below zero. The ternary inter-compound figure again shows that large negative values survive even for fairly dense structures, especially in systems such as `V-Zr-O`, `Sc-Zr-O`, `Nb-Ni-O`, `Ni-Ti-O`, `W-Y-O`, `Hf-Re-S`, `Pd-Zr-O`, `Mo-Zr-S`, and `V-Zr-S`.

The ternary within-chemsys histogram is broader than the binary one and shows strong populations near both `-1` and `+1`. That is a sign that many ternary chemsys have only a small number of available configurations, so some extreme correlation values are numerically sharp but chemically fragile. Even with that caution, the distribution still supports a broad negative packing association.

### Main takeaway

Packing efficiency has a clearer global relationship with `delta_E` than electronegativity does, but the relationship is not strong enough to explain the full chemistry. Dense systems can still be strongly negative, especially in the ternary space.

## 2. Maximum-delta_E point analysis

### Figures

![Binary maximum delta_E points](least-negative-delta-snapshots/plots/binary_max_delta_points.png)

![Ternary maximum delta_E points](least-negative-delta-snapshots/plots/ternary_max_delta_points.png)

### Interpretation

These plots retain only the numerically largest `delta_E` value for each category. Since the dataset is dominated by negative values, this corresponds to the least-negative or closest-to-zero point available for each binary metal or ternary metal pair.

The binary figure is tightly clustered near zero for most metals, showing that many binary systems have at least one strain point where the dip is weak. However, a few systems remain clearly negative even in that best-case view. The strongest examples are `Cu-O`, `Rh-S`, `Mo-O`, `Tc-O`, and `Zn-O`.

The ternary figure is visibly broader and more negative. Many ternary systems remain below about `-0.05`, and a smaller but important set stays below `-0.3` or even `-0.5`. Persistent examples include `Pd-Zr-O`, `Mo-Zr-S`, `V-Zr-S`, `Fe-Zr-S`, `Mo-Zr-O`, `Re-Zr-O`, `W-Zr-O`, `W-Zr-S`, and `Hf-V-S`.

### Main takeaway

Binary systems often have a near-zero best-case point. Ternary systems, in contrast, frequently remain appreciably negative even after selecting the most favorable available point. This is one of the clearest differences between the two datasets.

## 3. Grouped bar statistics

### Representative figures

![Binary oxygen min delta_E bars](aggregate-comparison-bars/plots/binary_oxygen_min_delta_bar.png)

![Binary sulfur min delta_E bars](aggregate-comparison-bars/plots/binary_sulfur_min_delta_bar.png)

![Ternary oxygen min delta_E bars](aggregate-comparison-bars/plots/ternary_oxygen_min_delta_bar.png)

![Ternary sulfur min delta_E bars](aggregate-comparison-bars/plots/ternary_sulfur_min_delta_bar.png)

![Binary oxygen mean delta_E bars](aggregate-comparison-bars/plots/binary_oxygen_mean_delta_bar.png)

### Interpretation

The `min delta_E` bar plots show the worst-case ranking most clearly.

For binary oxygen, `Ti-O` is the dominant outlier at about `-1.606`. The next tier includes `V-O`, `Y-O`, `Mn-O`, `Fe-O`, and `Cu-O`, but none reach the same severity. For binary sulfur, the deepest outlier is `Nb-S` at about `-1.084`, followed by `Zr-S`, `Ti-S`, and `V-S`.

The ternary oxygen `min` plot is led very clearly by `V-Zr-O` at about `-1.683`, which visually separates from the rest of the pair set. The next severe group includes `Sc-Zr-O`, `Nb-Ni-O`, `Mn-Ni-O`, `W-Y-O`, and `Ni-Ti-O`. In ternary sulfur, the deepest tail contains systems such as `Hf-Re-S`, `Mo-Zr-S`, `V-Zr-S`, `Fe-Zr-S`, and other transition-metal combinations involving Zr or Hf.

The mean and median bar plots compress much of the spread seen in the `min` panels. That means some systems are driven by deep excursions rather than by uniformly negative behavior at every sampled point. The binary oxygen mean plot illustrates this well: `Cu-O` remains the most negative because it has only one configuration, while `Ti-O`, `Y-O`, `Hf-O`, `Zn-O`, and `Sc-O` still stay noticeably below the rest.

The `max` bar families sit much closer to zero, which is consistent with the maximum-point analysis. Most groups have at least one comparatively mild point, but a small number remain negative even at their maximum retained value.

### Main takeaway

The grouped bar plots are the clearest way to see which systems dominate the negative tail. They show that ternary outliers are both more numerous and more persistent, and that the worst-case ordering is not always the same as the mean or median ordering.

## 4. Binary electronegativity trends

### Single-regression figures

![Binary oxygen min electronegativity plot](pauling-trend-checks/plots/binary_oxygen_min_electronegativity_vs_delta.png)

![Binary sulfur max electronegativity plot](pauling-trend-checks/plots/binary_sulfur_max_electronegativity_vs_delta.png)

### Multi-regression figure

![Binary electronegativity multi-regression grid](layered-pauling-regressions/plots/binary_clean_statistic_comparison_grid.png)

### Interpretation

The electronegativity figures show some signal, but not nearly as much structural organization as the packing plots. The clearest oxygen trend appears in `min delta_E`, where the regression is moderately positive. Lower-electronegativity metals such as `Y`, `Sc`, and especially `Ti` lie deeper in the negative region, whereas higher-electronegativity metals such as `Mo` and `Rh` are less negative on the `min` axis.

The clearest sulfur trend appears in `max delta_E`, where the regression is distinctly negative. In that panel, higher-electronegativity metals such as `Rh`, `Mo`, and `Pd` retain more negative maximum values, while several lower-electronegativity systems remain closer to zero.

The multi-regression grid is especially helpful because it shows that the strongest trends are concentrated in the extreme statistics rather than in the center of the distribution:

- `O min`: moderate positive trend
- `S min`: weak positive trend
- `O max`: weak negative trend
- `S max`: strongest negative trend
- `mean` and `median`: nearly flat for both ligands

This suggests that electronegativity is more relevant to the limiting behavior of the distribution than to the average response.

### Main takeaway

Electronegativity is a secondary descriptor in this dataset. It contributes ligand-dependent trends, especially in the extreme statistics, but it does not explain the full pattern of outliers.

## 5. Overall picture

Taken together, the figure families support a consistent picture:

- packing efficiency has a stronger global association with `delta_E` than electronegativity does
- severe negative outliers persist even at medium or high packing fraction
- ternary systems are more difficult than binary systems in both worst-case and best-case views
- Zr-containing ternaries appear repeatedly among the strongest negative outliers
- the `min` statistics carry the strongest chemical separation, while mean and median values compress much of the contrast

## 6. Short report summary

One concise way to describe the full set of results is:

> Packing efficiency shows a broad negative association with `delta_E`, but the relationship is not sufficient to explain the full behavior because strongly negative values remain even in relatively dense systems. Binary systems often recover to near-zero best-case values, whereas ternary systems retain more severe and more persistent negative dips. Electronegativity contributes weaker, ligand-dependent trends that are most visible in the extreme statistics rather than in the mean or median response.

## 7. Natural next steps

- trace the repeated Zr-containing ternary outliers back to their structural motifs
- compare packing fraction with volume-per-atom and coordination descriptors
- check whether the strongest ternary outliers are dominated by a small number of configuration families
- test robust regression for electronegativity so the `min` and `max` fits are less sensitive to outliers

