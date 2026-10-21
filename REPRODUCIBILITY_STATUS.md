# Computational reproducibility status

The repository accompanies **Gut microbial resource use can be low-dimensional yet functionally complex**.

## Checks executed

- **15** Python tests pass.
- **25** machine-readable expected-result tables parse correctly.
- Stored-results validation passes the main numerical checks for participation dimension, hCom2* thresholds, Ho support advantage, and rank-calibration thresholds.
- All **seven** data-driven figure sets rebuild from their machine-readable source tables into PNG and PDF. Figure 1 is provided as a conceptual graphic.
- Synthetic-input tests execute ordinary NMF, support-preserving strengths, species-identity nulls, and degree-preserving support rewires.
- The synthetic consumer–resource verification shows inverse-timescale tracking-error decay; exact tabulated Monte Carlo quantiles are preserved as result tables rather than represented as a verbatim random-draw replay.

## Biological source-model execution

The full hCom2* ecological prediction engine and Ho et al. MATLAB assembly prediction engine are distributed by their source-study authors. They are **not included** as newly authored Python implementations. Public sources are pinned and can be acquired with `scripts/fetch_sources.py`; the repository includes source-compatible matrix reconstruction and support-null generation, reference outputs, and the figure/claim mapping. Fully independent regeneration of all biological prediction results requires executing those source programs on the corresponding public experimental data.

### Interpretation boundaries

- hCom2* passage three is a fixed benchmark influenced by the source-study parameter choices, not a wholly prospectively held-out outcome.
- Ho species-label shuffling is the primary valid randomized comparison; degree-preserving rewiring is exploratory because some fitted models do not meet the feasibility criterion.
- The Ho source model constrains positive consumption rates to be equal within species.
- Initializing optimizers or running floating-point kernels with different library versions can affect exact numerical agreement.

No analysis results are claimed to have been computed from third-party source files in this build environment. These inputs are not redistributed.
