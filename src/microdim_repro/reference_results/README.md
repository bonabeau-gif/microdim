# Reference results

These CSV files are the numerical results printed in the final manuscript and supplementary information, plus code-preserved intermediate output from the sample-size calibration.

- `table_S4_cohort_benchmark.csv`: eight-cohort RRR/PCR/ridge results.
- `table_S5_niche_dimensions.csv`: abundance-based niche dimensionality and metabolite-space dimensionality.
- `table_S6_detection_thresholds.csv`: sample-size-dependent rank-10 detection thresholds.
- `table_S7_feature_cap.csv`: preprocessing-cap sensitivity.
- `table_S8_subspace_recovery.csv`: selected dimension fraction versus captured true subspace.
- `table_S9_micom_dimensions.csv`: MICOM model-implied consumption dimensions.
- `sample_matched_primary.csv`, `sample_and_cap_matched_sensitivity.csv`, `spectrum_sensitivity.csv`: exact original calibration outputs.
- `fast_resource_reported_summary.csv`: values stated in the manuscript from the original fast-resource run.
- `fast_resource_reconstructed_*`: output of the fully specified reconstructed numerical experiment distributed in this package.

These tables are intentionally kept separate from raw-data analysis functions. They are regression/audit targets, not hidden inputs to the analyses.
