# Ho et al. resource-sharing support analysis

## Public source model

Input archive: https://doi.org/10.5281/zenodo.7535703, including `PredictAssemblies_BHI.m`, `CR.m`, support/resource `.mat` arrays, fitted experimental resource levels, and the assembly outcome data. Checksum information is in `sources.lock.json`.

## Primary outcome

The primary external set consists of assembly identifiers **106–169**, inclusive (64 communities of 3–13 species). Model parameters and control variants are fitted only from permitted monoculture/spent-medium inputs, not the primary evaluation outcomes.

## Comparisons

- Experimentally inferred published support (benchmark).
- Species-specific resources only.
- Species-specific + pairwise resource sharing.
- Species-specific + all-but-one resources.
- Species-identity-shuffled support (100 computational realizations).
- Degree-preserving rewiring (100 realizations, 75 meeting resource-fit validity criteria; exploratory).

For each valid model compute the mean published prediction-error statistic across primary assemblies and the fraction of assemblies in which published support has lower error. The published support gives mean error **1.725205**, versus **2.383915** for identity shuffling, a **27.6%** improvement and **56/64** wins. These summaries are in `expected/ho_support_controls.csv`.

The published computational system constrains positive rates to be equal within each species. Consequently the support-map comparison tests support identity; it does not supply an independent empirical test of heterogeneous conditional strength compression.

## Commands

```bash
python scripts/fetch_sources.py --ho
python scripts/run_ho_validation.py
```

The validation script records checksums and the primary assembly IDs. Execute the provided MATLAB model to generate full assembly/control predictions; the Python archive checks the metadata and exposes the reference statistics but does not silently substitute a synthetic community model.

## Reproducible support-null generation

Export the published support array with species as rows and resource groups as columns, then run:

```bash
python scripts/run_support_nulls.py --support-csv /path/to/ho_support.csv --replicates 100
```

The script produces species-label permutations, degree-preserving double-edge swaps, and a manifest containing every random seed and accepted swap count. Each candidate support must independently undergo the source model's permitted parameter fitting and feasibility checks. Generated nulls alone are not predictions or valid community-model realizations.
