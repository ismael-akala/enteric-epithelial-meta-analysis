# ==============================================================================
# 01_load_and_filter_metadata.R
# 5-Cohort IBD Transcriptomic Meta-Analysis
# - Deterministic GSM-level master manifest verification (458 total GSMs: 368 included, 90 excluded)
# - Explicit platform verification (GPL570 / GPL6244) on getGEO() ExpressionSets
# - Per-cohort documented preprocessing & scale verification 
# - Probe-level audit counts before/after symbol mapping and limma::avereps collapsing
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(dplyr)
  library(tibble)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# Load verified GSM-level master manifest
master_manifest <- read.csv("verified_gsm_master_manifest.csv", stringsAsFactors = FALSE)

stopifnot(
  all(c("GSM", "GSE", "Status", "Final_Group", "Original_Label", "Source_Name", "Classification_Basis", "Exclusion_Reason") %in% colnames(master_manifest)),
  nrow(master_manifest) == 458L,
  length(unique(master_manifest$GSM)) == 458L,
  sum(master_manifest$Status == "Included") == 368L,
  sum(master_manifest$Status == "Excluded") == 90L
)

# Cohort specifications with explicit expected platform and documented GEO scale history
cohort_specs <- list(
  GSE16879 = list(
    platform = "GPL570",
    expected_ctrl = 12L, expected_ibd = 61L, expected_excl = 60L,
    scale_policy = "linear_to_log2",
    preprocessing_note = "Affymetrix HG-U133 Plus 2.0 (GPL570); GEO series matrix is on a linear intensity scale (MAS5/linear series matrix; 99th pct ~4581, max ~53971), so log2(pmax(x, 0) + 1) is applied prior to quantile normalization."
  ),
  GSE75214 = list(
    platform = "GPL6244",
    expected_ctrl = 22L, expected_ibd = 172L, expected_excl = 0L,
    scale_policy = "log2_already",
    preprocessing_note = "Affymetrix Human Gene 1.0 ST (GPL6244); GEO series matrix RMA log2-normalized in original series; no additional log transform or negative clipping applied."
  ),
  GSE38713 = list(
    platform = "GPL570",
    expected_ctrl = 13L, expected_ibd = 30L, expected_excl = 0L,
    scale_policy = "log2_already",
    preprocessing_note = "Affymetrix HG-U133 Plus 2.0 (GPL570); GEO series matrix RMA log2-normalized in original series; no additional log transform or negative clipping applied."
  ),
  GSE36807 = list(
    platform = "GPL570",
    expected_ctrl = 7L, expected_ibd = 28L, expected_excl = 0L,
    scale_policy = "log2_already",
    preprocessing_note = "Affymetrix HG-U133 Plus 2.0 (GPL570); GEO series matrix RMA log2-normalized in original series; no additional log transform or negative clipping applied."
  ),
  GSE4183 = list(
    platform = "GPL570",
    expected_ctrl = 8L, expected_ibd = 15L, expected_excl = 30L,
    scale_policy = "linear_to_log2",
    preprocessing_note = "Affymetrix HG-U133 Plus 2.0 (GPL570); GEO series matrix is on a linear intensity scale (99th pct ~11578, max ~70571), so log2(pmax(x, 0) + 1) is applied prior to quantile normalization."
  )
)

parse_hgnc_symbol <- function(annot_df, platform_id) {
  if (platform_id == "GPL570") {
    sym_col <- grep("^Gene Symbol$|^SYMBOL$", colnames(annot_df), ignore.case = TRUE, value = TRUE)[1]
    stopifnot(!is.na(sym_col))
    raw_sym <- trimws(as.character(annot_df[[sym_col]]))
    clean_sym <- sub("\\s*///.*$", "", raw_sym)
  } else if (platform_id == "GPL6244") {
    ga_col <- grep("gene_assignment|^gene symbol$|^symbol$", colnames(annot_df), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(ga_col)) {
      stop(sprintf("[GPL6244] Could not find gene_assignment or symbol column in annotation table. Available columns: %s", paste(colnames(annot_df), collapse = ", ")))
    }
    raw_ga <- as.character(annot_df[[ga_col]])
    if (grepl("^gene_assignment$", ga_col, ignore.case = TRUE)) {
      # Parse primary HGNC symbol from transcript block: RefSeq // SYMBOL // Description // ...
      clean_sym <- vapply(strsplit(raw_ga, "///"), function(blocks) {
        for (blk in blocks) {
          parts <- trimws(strsplit(blk, "//")[[1]])
          if (length(parts) >= 2L && nzchar(parts[2]) && parts[2] != "---") {
            return(parts[2])
          }
        }
        return("")
      }, character(1))
    } else {
      clean_sym <- sub("\\s*///.*$", "", trimws(raw_ga))
    }
  } else {
    stop("Unsupported platform: ", platform_id)
  }
  clean_sym[is.na(clean_sym) | clean_sym %in% c("", "---", "NA")] <- NA_character_
  return(clean_sym)
}

verify_expression_scale <- function(expr_mat, gse_id, policy) {
  q <- quantile(expr_mat, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  n_neg <- sum(expr_mat < 0, na.rm = TRUE)
  # Check intensity range consistency for log2-scale GEO series.
  if (policy == "log2_already") {
    if (q["99%"] > 25 || q["100%"] > 30) {
      stop(sprintf("[%s] Expected log2-scale GEO matrix, but 99th percentile=%.2f, max=%.2f. Inspect series preprocessing.", gse_id, q["99%"], q["100%"]))
    }
    return(list(mat = expr_mat, transformed = FALSE, q0 = q["0%"], q50 = q["50%"], q99 = q["99%"], q100 = q["100%"], n_neg = n_neg))
  } else if (policy == "linear_to_log2") {
    if (n_neg > 0L) {
      message(sprintf("[%s] Linear series contains %d negative values (min=%.4f); applying log2(pmax(x, 0) + 1).", gse_id, n_neg, q["0%"]))
    }
    mat_log2 <- log2(pmax(expr_mat, 0) + 1)
    return(list(mat = mat_log2, transformed = TRUE, q0 = q["0%"], q50 = q["50%"], q99 = q["99%"], q100 = q["100%"], n_neg = n_neg))
  }
  stop("Unknown scale policy: ", policy)
}

cohort_data_list <- list()
preprocessing_audit <- list()
probe_audit <- list()

for (gse_id in names(cohort_specs)) {
  spec <- cohort_specs[[gse_id]]
  message(sprintf("[%s] Retrieving GEO series and verifying platform %s...", gse_id, spec$platform))
  
  gset_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
  if (length(gset_list) == 0L) stop("No ExpressionSet returned for ", gse_id)
  
  # Explicitly select and verify the ExpressionSet matching the expected platform
  platforms_found <- vapply(gset_list, annotation, character(1))
  match_idx <- which(platforms_found == spec$platform)
  if (length(match_idx) != 1L) {
    stop(sprintf("[%s] Expected exactly 1 ExpressionSet on %s, found: %s", gse_id, spec$platform, paste(platforms_found, collapse = ", ")))
  }
  eset <- gset_list[[match_idx]]
  
  # Subset strictly against verified GSM manifest
  sub_manifest <- master_manifest %>% filter(GSE == gse_id)
  inc_manifest <- sub_manifest %>% filter(Status == "Included")
  exc_manifest <- sub_manifest %>% filter(Status == "Excluded")
  
  stopifnot(
    sum(inc_manifest$Final_Group == "Control") == spec$expected_ctrl,
    sum(inc_manifest$Final_Group == "IBD") == spec$expected_ibd,
    nrow(exc_manifest) == spec$expected_excl,
    all(sub_manifest$GSM %in% sampleNames(eset))
  )
  
  eset_inc <- eset[, inc_manifest$GSM]
  stopifnot(identical(sampleNames(eset_inc), inc_manifest$GSM))
  
  expr_raw <- exprs(eset_inc)
  n_probes_raw <- nrow(expr_raw)
  
  # Per-cohort scale verification.
  scale_res <- verify_expression_scale(expr_raw, gse_id, spec$scale_policy)
  expr_norm <- limma::normalizeBetweenArrays(scale_res$mat, method = "quantile")
  
  # Map probes to primary HGNC gene symbols and collapse replicates
  fdat <- fData(eset_inc)
  symbols <- parse_hgnc_symbol(fdat, spec$platform)
  valid_idx <- which(!is.na(symbols))
  n_probes_mapped <- length(valid_idx)
  
  expr_mapped <- expr_norm[valid_idx, , drop = FALSE]
  rownames(expr_mapped) <- symbols[valid_idx]
  expr_collapsed <- limma::avereps(expr_mapped, ID = rownames(expr_mapped))
  n_genes_collapsed <- nrow(expr_collapsed)
  
  cohort_data_list[[gse_id]] <- list(
    expr = expr_collapsed,
    pheno = inc_manifest,
    platform = spec$platform
  )
  
  preprocessing_audit[[gse_id]] <- data.frame(
    GSE = gse_id,
    Platform = spec$platform,
    Total_Series_GSMs = ncol(eset),
    Retained_GSMs = nrow(inc_manifest),
    Control_n = sum(inc_manifest$Final_Group == "Control"),
    IBD_n = sum(inc_manifest$Final_Group == "IBD"),
    Excluded_GSMs = nrow(exc_manifest),
    Scale_Policy = spec$scale_policy,
    Log2_Transformed_In_Pipeline = scale_res$transformed,
    Min_Intensity = round(scale_res$q0, 4),
    Median_Intensity = round(scale_res$q50, 4),
    P99_Intensity = round(scale_res$q99, 4),
    Max_Intensity = round(scale_res$q100, 4),
    Negative_Values_Count = scale_res$n_neg,
    Preprocessing_Documentation = spec$preprocessing_note,
    stringsAsFactors = FALSE
  )
  
  probe_audit[[gse_id]] <- data.frame(
    GSE = gse_id,
    Platform = spec$platform,
    Raw_Probes = n_probes_raw,
    Mapped_Probes = n_probes_mapped,
    Unmapped_Probes = n_probes_raw - n_probes_mapped,
    Collapsed_Unique_Genes = n_genes_collapsed,
    stringsAsFactors = FALSE
  )
}

# Strict 5-cohort gene symbol intersection
common_genes <- sort(Reduce(intersect, lapply(cohort_data_list, function(x) rownames(x$expr))))
message(sprintf("Strict 5-cohort intersection yielded %d common HGNC gene symbols.", length(common_genes)))

for (gse_id in names(cohort_data_list)) {
  cohort_data_list[[gse_id]]$expr <- cohort_data_list[[gse_id]]$expr[common_genes, , drop = FALSE]
  stopifnot(identical(rownames(cohort_data_list[[gse_id]]$expr), common_genes))
}

# Export auditable manifests and preprocessing logs
write.csv(master_manifest %>% filter(Status == "Included"), "results/reconciled_included_manifest.csv", row.names = FALSE)
write.csv(master_manifest %>% filter(Status == "Excluded"), "results/reconciled_exclusion_log.csv", row.names = FALSE)
write.csv(bind_rows(preprocessing_audit), "results/cohort_preprocessing_audit.csv", row.names = FALSE)
write.csv(bind_rows(probe_audit) %>% mutate(Intersected_Common_Genes = length(common_genes)), "results/probe_collapsing_audit.csv", row.names = FALSE)

saveRDS(cohort_data_list, "results/harmonized_5cohorts.rds")
saveRDS(common_genes, "results/common_genes.rds")
capture.output(sessionInfo(), file = "results/sessionInfo_01.txt")
