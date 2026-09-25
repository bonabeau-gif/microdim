# Reproducibility notes and provenance

## 1. Why the package contains frozen reference tables

The manuscript was assembled from several analysis stages. The later sample-size calibration and MICOM extension retained executable scripts and/or fully public raw model output. Earlier observational and numerical stages retained complete scientific descriptions and final numerical summaries, but not every incidental software choice. Rather than invent a history, this package stores the published values as reference targets and labels reimplemented analyses as such.

This is preferable to hard-coding the published answers into the analysis routines: the routines operate on raw data, while the frozen tables provide an audit target.

## 2. Exact components

### Rank calibration

The files in `docs/original_rank_calibration_script.py` and `docs/original_rank_validation_script.py` are preserved copies of the scripts available when the final calibration was performed. The package module is a refactor of the same formulas and seeds.

### MICOM

The source exchange table is public and contains the exact flux values used for the manuscript extension. The package computes its metrics directly from those rows.

## 3. Reimplemented components

### Cohort benchmark

The manuscript uniquely fixes the outer-training-only preprocessing logic, nested group-aware validation, model classes, and endpoint. It does not uniquely fix the historical pseudocount constant, random fold seed, or all tuning grids. The package makes these explicit. For strict comparison, use `reference_results/table_S4_cohort_benchmark.csv`.

### Effective niche dimensionality

The published estimator fixes the likelihood, gradients, dimensions, reconstruction metric, and stopping criterion. Initialization and optimizer schedule can change finite-sample fits. The package uses deterministic initialization and backtracking. The manuscript values remain in `table_S5_niche_dimensions.csv`.

### Fast-resource numerical check

The exact historical random parameter file was not retained. The reconstruction uses stated dimensions, factorized positive uptake, 30 deterministic seeds, the manuscript eta grid, quasi-steady initialization, 101 common time points, and tight ODE tolerances. It independently verifies the theoretical signatures. Manuscript-reported summary values are preserved separately.

## 4. What counts as reproduction

There are three useful levels:

1. **Figure/table regeneration** from the published numerical outputs (`microdim-repro figures`).
2. **Algorithmic reproduction** from public raw inputs using fully documented current code.
3. **Historical bitwise reproduction**, possible only where original scripts/seeds/input files survived.

The package explicitly identifies which level applies to each result.
