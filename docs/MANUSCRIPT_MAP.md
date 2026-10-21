# Manuscript claim-to-code and claim-to-table map

| Manuscript analysis | Code | Machine-readable result |
| --- | --- | --- |
| Fig. 2: Cohort held-out prediction | `src/microdim/cohort.py` | `expected/cohort_benchmark.csv` |
| Supplementary Table 2 | `scripts/run_cohort_benchmark.py` | `expected/cohort_benchmark.csv` |
| Supplementary Fig. 1: niche dimensions | `src/microdim/niche.py` | `expected/niche_dimensions.csv` |
| Supplementary Table 3 | same | `expected/niche_dimensions.csv` |
| Fig. 3, Supplementary Tables 4–7 | `src/microdim/calibration.py` | `expected/rank_calibration.csv`, `rank_thresholds.csv`, `subspace_recovery.csv` |
| Reduced consumer–resource validity (Supplementary Table 1) | `src/microdim/cr_theory.py` | `expected/cr_timescale.csv` |
| hCom2* baseline (Supplementary Table 8) | Published model kernel | `expected/hcom2_baseline.csv` |
| Consumption geometry (Supplementary Table 9) | `src/microdim/metrics.py` | `expected/hcom2_geometry.csv` |
| Fig. 4, Table 1, Supplementary Table 10 | `src/microdim/support_strength.py`, source prediction kernel | `expected/hcom2_nmf_threshold.csv`, `hcom2_nmf_selected_highrank.csv` |
| Fixed support, Table 1, Supplementary Table 11 | `src/microdim/support_strength.py`, source prediction kernel | `expected/hcom2_support_strength.csv`, `hcom2_support_selected_full.csv` |
| Supplementary Fig. 3: weak-edge removal | Source model kernel | `expected/hcom2_weak_edge_ablation.csv` |
| Fig. 5, Table 2 | Ho source prediction kernel | `expected/ho_support_controls.csv`, `ho_support_controls_full.csv` |
| Other numerical robustness diagnostics | Relevant workstreams | `expected/hcom2_*`, `expected/ho_*` |

Every expected table is a frozen quantitative snapshot, not a substitute for running the relevant third-party source computation.
