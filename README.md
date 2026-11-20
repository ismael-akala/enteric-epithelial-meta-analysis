# Enteric Epithelial Transcriptomic Meta-Analysis (5-Cohort IBD)

Reproducible R pipeline and meta-analytic workflow for identifying consensus transcriptomic signatures across five public Inflammatory Bowel Disease (IBD) cohorts (`GSE16879`, `GSE75214`, `GSE38713`, `GSE36807`, and `GSE4183`; 368 included samples: 306 IBD and 62 controls, with 90 excluded).

## Repository Structure

- `scripts/`
  - `01_load_and_filter_metadata.R` — Verifies GSM manifest counts, checks GEO platform/intensity scales, performs quantile normalization, maps probes to HGNC symbols, collapses replicates, and extracts the 5-cohort gene intersection.
  - `02_limma_differential_expression.R` — Runs two-arm `limma` models (`IBD` vs. `Control`) per cohort and extracts moderated standard errors.
  - `03_metafor_rema.R` — Fits DerSimonian–Laird random-effects models with Knapp–Hartung adjustment (`df = 4`) across intersected genes (`FDR < 0.05`, `|Pooled_log2FC| >= 0.50` for consensus DEGs).
  - `04_clusterProfiler_enrichment.R` — Performs directional GO Biological Process (`GO:BP`) enrichment with Entrez mapping audit logs.
- `data/` — Verified sample manifest (`verified_gsm_master_manifest.csv`).
- `results/` — Cohort preprocessing/probe/Entrez audit logs, meta-analysis tables, and GO:BP enrichment outputs.

## Archival & Citation

