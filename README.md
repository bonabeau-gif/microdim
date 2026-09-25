# Reproducibility package: predictive, niche, and mechanistic dimensionality

This package reproduces the computational analyses for the manuscript
**“What dimension have we estimated? Dimensionality in microbiome-metabolome models.”**

It covers:

1. the eight-cohort reduced-rank regression (RRR), principal-components regression (PCR), and ridge benchmark;
2. the abundance-based effective niche dimensionality analysis following Srinivasan, Plata & Dixit (2026);
3. the sample-size-matched rank-selection calibration and its sensitivities;
4. the consumer-resource factorization and fast-resource numerical check;
5. the MICOM model-implied consumption-dimensionality analysis;
6. all main and supplementary figures.

## Important reproducibility status

Not all historical analyses have the same provenance. This package makes that explicit.

**Exact code-preserved analyses.** The sample-size-matched rank calibration is incorporated from the original script used for the manuscript, including seeds, grids, and the ridge penalty. The MICOM analysis is fully reproducible from a public exchange-flux file and the package downloads that file directly.

**Public-data re-executable analyses.** The eight-cohort benchmark and effective-niche analysis can be rerun from the public Borenstein Lab curated collection. The manuscript preserved the preprocessing rules and model definitions, but not every legacy random seed, optimization step size, or pseudocount implementation detail. The Python implementation here fixes those choices explicitly and supplies the manuscript tables as regression targets in `reference_results/`. Small numerical differences in a from-scratch rerun are therefore possible and should not be silently interpreted as new biological results.

**Reconstructed numerical check.** The manuscript preserved the dimensions, model equations, timescale grid, number of seeds, and reported summary statistics for the fast-resource simulation, but the exact original random parameter ranges and seed file were not retained. `consumer_resource.py` implements a deterministic reconstruction of the stated experiment and separately preserves the manuscript-reported numbers in `reference_results/fast_resource_reported_summary.csv`. The reconstruction verifies the same rank, sign, reciprocity, positivity, and approximately inverse-timescale error behavior without pretending to recover lost historical random draws.

This separation is deliberate: reference-result regeneration and independent re-execution are both provided, but they are not conflated.

## Installation

With conda:

```bash
conda env create -f environment.yml
conda activate microdim-repro
```

Or with pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## One-command workflow

```bash
microdim-repro all --data-dir data --results-dir results
```

This downloads the required public tables and runs all analyses. The effective-niche fit is the slowest step because it fits dimensions 1–5 across 100 subsamples for eight cohorts.

For a fast installation check:

```bash
pytest -q
microdim-repro figures --outdir results/reference_figures
```

## Individual analyses

```bash
microdim-repro download --data-dir data
microdim-repro cohorts --data-dir data --outdir results/cohorts
microdim-repro niche --data-dir data --outdir results/niche --subsamples 100
microdim-repro niche-sensitivity --data-dir data --outdir results/niche_sensitivity --subsamples 100
microdim-repro calibration --outdir results/calibration
microdim-repro consumer-resource --outdir results/consumer_resource
microdim-repro micom --data-dir data --outdir results/micom
microdim-repro figures --outdir results/reference_figures
```

`figures` regenerates manuscript-style figures from the frozen manuscript result tables. It does not require downloading the large cohort tables. Plot titles are deliberately omitted because the manuscript captions carry the explanations.

## Public data downloaded

For each of the eight paired cohorts, only three processed files are downloaded from
`borenstein-lab/microbiome-metabolome-curated-data`:

- `genera.tsv`
- `mtb.tsv`
- `metadata.tsv`

The selected cohorts are:

- `ERAWIJANTARI_GASTRIC_CANCER_2020`
- `FRANZOSA_IBD_2019`
- `HE_INFANTS_MFGM_2019`
- `JACOBS_IBD_FAMILIES_2016`
- `KANG_AUTISM_2017`
- `WANDRO_PRETERMS_2018`
- `WANG_ESRD_2020`
- `YACHIDA_CRC_2019`

The MICOM analysis downloads:

- `urtecholabucsd/FMT/Figure1/input/exchanges_d0_high_fiber.csv`
- `urtecholabucsd/FMT/Figure1/processed_data/vendor_diversity.csv`

Third-party data retain their original licenses and terms. They are not redistributed in this ZIP.

## Analysis-specific notes

### Eight-cohort predictive benchmark

Every outer training fold learns its own preprocessing. Taxa must be present in at least 25% of training samples and metabolites observed in at least 75%. Both blocks are capped at one quarter of outer-training sample size. Taxa are zero-replaced, CLR transformed, and standardized. Metabolites are median-imputed, log-transformed, and standardized. Group-aware folds use family for Jacobs and subject otherwise.

The historical manuscript did not retain the exact pseudocount constant or every model-selection grid. This package uses one half of the minimum positive training value as a fold-specific pseudocount, inner-CV ridge alphas from `1e-3` to `1e3`, PCR dimensions 0–12, and RRR ranks 0–12 with ridge stabilization 10. These choices are documented in code and can be changed.

### Effective niche dimensionality

The implementation follows the published model

`x_so = softmax_o(sum_k z_sk theta_ko)`

with multinomial likelihood, the published gradients, dimensions 1–5, and the published relative-gradient stopping criterion. Each cohort is filtered at 0.1% mean relative abundance and subsampled 100 times to 50 samples, or 70% of the cohort if fewer than 50 are available. `eta_D` is the exponential decay scale of mean KL reconstruction error versus latent dimension.

The publication specifies the objective and stopping rule but not a unique random initialization and step-size schedule. The deterministic implementation therefore serves as a transparent reimplementation, while `reference_results/table_S5_niche_dimensions.csv` stores the manuscript values.

### Rank-selection calibration

The code is a package version of the original manuscript script. It uses 100 seeds per condition, true ranks 2 and 10, `p=q=24`, the cohort-matched rank-selection sample sizes, a fixed ridge penalty of 10, and population R² grid `0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50`. It also reproduces feature-cap and singular-spectrum sensitivities and Table S8 subspace recovery.

### Consumer-resource model

The module implements

- logistic self-renewing resources;
- bilinear uptake;
- the quasi-steady resource manifold;
- the reduced Lotka–Volterra interaction matrix `A = -c diag(v/s) c.T`;
- the feature factorization `c = Phi_n C Phi_R.T`;
- numerical checks of factorization, rank, negative semidefiniteness, and finite-timescale tracking.

### MICOM proxy

For each community the script removes the `medium` rows, keeps `direction == "import"`, and constructs a taxon × metabolite matrix. It reports numerical rank, participation-ratio dimension, entropy effective rank, d90, d95, binary uptake dimension, and abundance-weighted dimension. The primary cutoff is `1e-6`, with sensitivity through `1e-3`.

## Reference results

`reference_results/` contains the numerical values printed in the final manuscript/supplement, plus the original calibration output files. These let users distinguish a code/version difference from a manuscript transcription error.

## Tests

The unit tests check:

- exact reduced-rank construction;
- the participation-ratio formula;
- the consumer-resource factorization/rank/sign constraints;
- calibration determinism for a small fixed condition;
- figure generation from reference results.

## Computational requirements

The reference figures and reconstructed consumer-resource check run on a laptop. The full rank calibration can take several minutes. The full effective-niche rerun is substantially slower because of the nested 8 × 100 × 5 optimization structure. Parallelization is intentionally not hidden inside the estimator; users can split cohorts if desired.
