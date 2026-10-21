# hCom2* consumer-resource prediction

## Input and kernel

The 63-strain, 98-resource input consumption array is built from the deposited Mahajan–Maslov source data. The source model provides resource concentration fitting, production normalization, and serial dilution; use the git commit pinned in `sources.lock.json` for source-level execution. The figures and numerical outcomes are recorded in `expected/hcom2_*.csv`.

## Analysis definitions

1. Calculate singular values of the nonnegative consumption matrix, participation dimension, entropy dimension, energy-retention dimensions, and numerical rank.
2. Fit ordinary NMF minimizing Frobenius reconstruction error at each tested rank with 20 deterministic initializations. The selected seed minimizes **matrix reconstruction error**, not passage-three prediction error.
3. Fit support-preserving strengths with a fixed measured binary mask. On positive entries the model is `sigmoid(a_i + b_j + u_i @ v_j)`; regularization and optimizer details are in `src/microdim/support_strength.py`.
4. Run a new source-model resource fit for every candidate reconstructed consumption matrix. Predict serial passages using the same kinetic program, source inputs and evaluation set.
5. Measure passage-three Pearson correlation and base-10 log-abundance RMSE. Use the monoculture null to compute retained improvement `(RMSE_null - RMSE_k)/(RMSE_null - RMSE_full)`.
6. Call a candidate sufficient when retained improvement is at least 0.95, Pearson degradation is no greater than 0.03, and the resource/prediction solver converges. Record seed-wise pass fractions.
7. For support robustness, compute weak-edge ablation, energy retention and random edge controls without changing the training and prediction protocol.

Selected-seed sufficiency requires ranks 52 (NMF) and 6 (fixed-support strengths); corresponding seed thresholds are stored in `expected/`.

## Commands

```bash
python scripts/fetch_sources.py --hcom2
python scripts/run_hcom2_analysis.py --consumption-csv /path/to/63x98_matrix.csv
python scripts/rebuild_figures_from_frozen.py
```

`run_hcom2_analysis.py` computes the measurable matrix-geometry diagnostics from an exported consumption matrix. Full passage predictions must be executed with the pinned published serial-dilution model. This repository does not replace its resource fit/prediction code.

**Interpretation:** The dimension six refers only to fitted *positive strengths conditional on the fixed support*. It is not the rank or complete complexity of the full resource-consumption matrix. Passage three is a benchmark used in the source modeling study's parameterization and should not be described as an untouched test.

## Independent reconstruction grid

```bash
python scripts/run_compression.py --consumption-csv /path/to/consumption.csv --ordinary-ranks 1,5,20,52,60 --support-dimensions 0,3,6,12 --seeds 0,1,2
```

This saves the complete numerical fitting grid and the reconstruction-selected consumption matrix at every requested dimension under `results/compression/`. The selection criterion is reconstruction loss alone. To obtain biological passage-three metrics, pass every reconstructed matrix through the published source model's independently fitted resource and serial-dilution prediction steps.
