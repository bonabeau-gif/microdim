# Execution and reproducibility coverage

All numerical claims should be read together with their analysis objects, evaluation criteria, and benchmark limitations. A stored numerical result is not evidence that its underlying biological simulation has been re-executed in the current environment.

| Workstream | Directly executable in this repository | Required outside inputs or programs |
| --- | --- | --- |
| Cohort RRR/PCR/ridge | `run_cohort_benchmark.py`, `src/microdim/cohort.py` | Paired processed cohort tables from pinned source repository |
| Latent niche dimensions | `run_niche_dimension.py`, `src/microdim/niche.py` | Same cohort input tables |
| Rank-detection calibration | `run_rank_calibration.py`, `src/microdim/calibration.py` | None; full 100-seed grid is compute intensive |
| Consumer–resource theory | `run_consumer_resource_verification.py`, `src/microdim/cr_theory.py` | None; reference random-system verification is not an exact draw replay of the tabulated quantiles |
| hCom2* matrix geometry | `run_hcom2_analysis.py` | Consumption matrix exported from published source |
| hCom2* passage predictions and refitting | Reference tables and core compression code | Mahajan–Maslov source-model kinetics and all biological inputs |
| Ho independent support control | Reference tables and fixed assembly-ID inventory | Published Ho MATLAB system and public inputs |
| Main figures 2–5 and SI figures 1–3 | `rebuild_figures_from_frozen.py` | None |

## Run without external sources

```bash
python scripts/run_all.py
python scripts/verify_tables.py
python scripts/verify_figures.py
```

## Run the Python analyses with public sources

```bash
python scripts/fetch_sources.py --cohorts --hcom2 --ho
python scripts/run_cohort_benchmark.py
python scripts/run_niche_dimension.py
python scripts/run_rank_calibration.py
python scripts/run_consumer_resource_verification.py
```

Use `docs/HCOM2.md` and `docs/HO.md` for model-specific source execution. Computations depending on external model code are attributed to those model programs; the machine-readable recorded outcomes permit quantitative manuscript auditing without rerunning them. Numerical random seeds are in `config/seeds.json`.

## Checks

- `pytest -q`: implemented algebra, regression, factorization, and archive assertions.
- `validate_frozen.py`: core scientific result thresholds.
- `verify_tables.py`: confirm all machine-readable data tables parse.
- `rebuild_figures_from_frozen.py`: regenerate numerical figures from their exact tabular inputs.
- `verify_figures.py`: confirm seven data-driven figure pairs are produced.

Do not interpret agreement with the `expected/` tables alone as an independent kinetic-model rerun. The deposited hCom2* test set is not an untouched prospective holdout; the Ho model imposes within-species positive-rate equality, and only the valid support-control parameterizations contribute to its exploratory randomized-family statistics.
