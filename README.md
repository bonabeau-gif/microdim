# Gut microbial resource use can be low-dimensional yet functionally complex

Research code and numerical reproducibility record for the manuscript by **Eric Bonabeau** (School of Complex Adaptive Systems, Arizona State University).

The analysis distinguishes spectral concentration, predictive rank, trait-mediated interaction rank, fixed support topology, and predictive sufficiency in microbial consumer–resource systems. It includes eight paired human microbiome–metabolome cohorts, a trait-factorized consumer–resource model, a 63-strain hCom2* serial-dilution system, and 15-species resource-sharing experiments from Ho et al.

## Get started

Python 3.11 is recommended (a Conda definition and a Dockerfile are included).

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
python scripts/validate_frozen.py
python scripts/verify_tables.py
python scripts/rebuild_figures_from_frozen.py
python scripts/verify_figures.py
```

These commands require no external downloads. They validate the manuscript's archived numerical results, exercise the analysis modules, and recreate Figures 2–5 and Supplementary Figures S1–S3 in `figures/rebuilt/`. Figure 1 is a conceptual diagram provided in `figures/reference/`.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/microdim/` | Quantitative metrics, regression, niche, calibration, consumer–resource and factorization code |
| `scripts/` | Testable experiment drivers, figure construction, numerical verification, source download |
| `expected/` | 25 numerical tables supporting the paper and supplement |
| `figures/reference/` | Submitted figure artwork |
| `figures/rebuilt/` | Outputs of deterministic figure-generation scripts |
| `config/` | Data mapping, optimizer seeds and tolerances |
| `sources.lock.json` | Commit-pinned public input/model repositories and checksums |
| `docs/` | Analyses, data access, dependencies and reproducibility scope |
| `.github/workflows/checks.yml` | GitHub Actions validation on Python 3.11 |

## Scientific analyses

**Predictive regression:** group-aware nested evaluation of RRR, PCR and ridge from paired cohorts (`scripts/run_cohort_benchmark.py`, `src/microdim/cohort.py`).

**Niche dimension:** multinomial latent-space fitting and reconstruction-decay estimation (`scripts/run_niche_dimension.py`, `src/microdim/niche.py`).

**Detection calibration:** rank-10 synthetic regression models, matched sample sizes, fixed seeds, and cross-validated rank selection (`scripts/run_rank_calibration.py`).

**Consumer–resource reduction:** resource quasi-steady reduction, positive-semidefinite interaction structure, and numerical timescale checks (`scripts/run_consumer_resource_verification.py`).

**hCom2*:** matrix geometry, nonnegative low-rank approximation and support-conditioned strengths (`scripts/run_hcom2_analysis.py`, `src/microdim/support_strength.py`). Resource refitting and serial-dilution prediction use the publicly deposited Mahajan–Maslov kinetic code.

**Ho et al.:** checksummed external model data, fixed 64-assembly outcome set, and comparison with alternative resource-sharing support structures. The archived MATLAB computational model supplies the prediction kernel.

See [methods and coverage](docs/EXECUTION.md), [claim-to-artifact map](docs/MANUSCRIPT_MAP.md), and [reproducibility status](REPRODUCIBILITY_STATUS.md).

## Running analyses with public biological data

```bash
python scripts/fetch_sources.py --cohorts --hcom2 --ho
python scripts/run_cohort_benchmark.py
python scripts/run_niche_dimension.py
python scripts/run_rank_calibration.py --smoke
python scripts/run_consumer_resource_verification.py
python scripts/run_hcom2_analysis.py --consumption-csv /path/to/consumption.csv
python scripts/run_ho_validation.py
```

Public data are placed under `data/third_party/` and are not included in this repository. The complete rank-calibration run is computationally intensive; omit `--smoke` to run 100 independent datasets per condition. Resource fits and community predictions require executing the external model programs described in `docs/HCOM2.md` and `docs/HO.md`. The Python entry points do **not** claim to independently execute those two kinetic systems.

## Reproducibility scope

The `expected/` tables are the machine-readable numerical record for the manuscript and provide an immediate route to check stated values and figures. Analytical and synthetic computations can be executed directly with the included Python modules. Some ecological evaluations require independently published source-model kernels and associated public input data. Where exact results depend on those kernels, the repository supplies pinned source identities, numerical tables, methods, and clear computational boundaries rather than substituting an unrelated model.

The passage-three hCom2* outcome informed model choices in the source study; it is a fixed computational benchmark, not a wholly untouched prospective test. Ho's primary support comparison uses the 64 preselected mixed assemblies. Its within-species equality of positive resource coefficients is imposed by the source model and does not independently establish heterogeneity-free biology.

## Manuscript result highlights

- hCom2* consumption-matrix participation dimension: **3.27**.
- Smallest tested functionally sufficient ordinary NMF rank: **52**, reconstruction-selected seed; support-preserving latent-strength dimension: **6**, with the binary support supplied separately.
- Seed-agreement thresholds: ordinary NMF **52 / 57 / 60** and fixed-support strengths **6 / 6 / 12** (selected / majority / all).
- Ho resource-sharing identity control: **27.6%** relative reduction in mean primary error for the published support; **56/64** assemblies favor it.

These are benchmark-specific operational quantities, not estimates of a universal niche count.

## Citation and licensing

Cite the accompanying manuscript and the original biological source studies as appropriate (see `CITATION.cff` and `docs/DATA_SOURCES.md`). Original repository Python code is licensed under MIT. Data and programs obtained from external repositories are governed by their respective licenses.

## Low-rank model grids and support nulls

The repository also includes `scripts/run_compression.py` to generate all reconstruction-selected ordinary-NMF and support-preserving candidate matrices, and `scripts/run_support_nulls.py` to generate reproducible species-identity and degree-preserving support randomizations. These outputs are inputs to the public ecological prediction kernels, not substitutes for them. See `docs/HCOM2.md` and `docs/HO.md`.
