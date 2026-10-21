# Public inputs

The source download and exact revision metadata are in `sources.lock.json`.

- **Paired cohorts:** Muller, Y. M. Algavi & Borenstein, *npj Biofilms and Microbiomes* (2022). Processed, paired taxonomic and metabolite tables from `borenstein-lab/microbiome-metabolome-curated-data`. Eight cohorts are enumerated in `config/cohorts.json`.
- **hCom2*:** Mahajan & Maslov, *PLOS Computational Biology* (2025), DOI 10.1371/journal.pcbi.1013222. Published Python code and deposited consumption/production data in `Tarun-Mahajan/crn_gut_microbiome`.
- **Experimental input sources:** Han et al., *Nature* (2021), DOI 10.1038/s41586-021-03707-9; Jin et al., *Nature Communications* (2023), DOI 10.1038/s41467-023-39121-0. Related repositories are pinned in `sources.lock.json`.
- **Ho model:** Ho et al., *Nature Microbiology* (2024), DOI 10.1038/s41564-024-01625-w. Deposited code and data: https://doi.org/10.5281/zenodo.7535703. Downloaded archive files with known checksums are checked by the fetcher.

Biological input files and external model source code are fetched into `data/third_party/`. They are excluded from source control; consult the original repositories for their data use and attribution terms. Download scripts require an internet connection.
