# Data sources

## Paired microbiome–metabolome cohorts

Muller E, Algavi YM, Borenstein E. 2022. The gut microbiome-metabolome dataset collection: a curated resource for integrative meta-analysis. npj Biofilms and Microbiomes 8:79.

Repository: `borenstein-lab/microbiome-metabolome-curated-data`

The package downloads processed genus tables, metabolite tables, and metadata only. It does not redistribute the underlying studies.

## MICOM flux proxy

Zhang Z, Holton M, Ferrer DM, Tripp AD, Richter A, Dixit PD, Urtecho G. 2026. Metagenome-scale Modeling to Assess Microbiome Metabolic Complementarity for Precision Microbiota Transplantation Therapies. bioRxiv. DOI: 10.64898/2026.05.15.725570.

Repository: `urtecholabucsd/FMT`

Used files:
- `Figure1/input/exchanges_d0_high_fiber.csv`
- `Figure1/processed_data/vendor_diversity.csv`

## MICOM software reference

Diener C, Gibbons SM, Resendis-Antonio O. 2020. MICOM: Metagenome-Scale Modeling To Infer Metabolic Interactions in the Gut Microbiota. mSystems 5:e00606-19.

## Pinned repository snapshots

The downloader uses immutable commit identifiers:

- Borenstein curated collection: `89a519d8c832008fbc6e650453e83e2f04858d02`
- FMT/MICOM exchange analysis: `f571e19221b18a8a2517d0a77fc30652ed6326b8`

This prevents a future repository update from silently changing the input data used by the reproduction workflow.
