# ==============================================================================
# 02_limma_differential_expression.R
# Within-cohort two-arm limma models (IBD vs Control) with direct moderated SE
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})

cohort_data_list <- readRDS("results/harmonized_5cohorts.rds")
common_genes <- readRDS("results/common_genes.rds")

de_results_list <- list()

for (gse_id in names(cohort_data_list)) {
  expr_mat <- cohort_data_list[[gse_id]]$expr
  pheno_df <- cohort_data_list[[gse_id]]$pheno

  stopifnot(
    identical(colnames(expr_mat), pheno_df$GSM),
    identical(rownames(expr_mat), common_genes),
    sum(pheno_df$Final_Group == "Control") >= 2L,
    sum(pheno_df$Final_Group == "IBD") >= 2L
  )

  group <- factor(pheno_df$Final_Group, levels = c("Control", "IBD"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("Control", "IBD")

  fit <- lmFit(expr_mat, design)
  contrast_mat <- makeContrasts(IBD_vs_Control = IBD - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)

  # Direct moderated standard error extraction
  coef_idx <- which(colnames(fit2$coefficients) == "IBD_vs_Control")
  stopifnot(length(coef_idx) == 1L)

  logFC <- as.numeric(fit2$coefficients[, coef_idx])
  SE <- as.numeric(fit2$stdev.unscaled[, coef_idx] * sqrt(fit2$s2.post))
  t_stat <- as.numeric(fit2$t[, coef_idx])
  p_val <- as.numeric(fit2$p.value[, coef_idx])

  stopifnot(
    all(is.finite(logFC)),
    all(is.finite(SE) & SE > 0)
  )

  res_df <- data.frame(
    GeneSymbol = rownames(expr_mat),
    GSE = gse_id,
    logFC = logFC,
    SE = SE,
    t = t_stat,
    P.Value = p_val,
    adj.P.Val = p.adjust(p_val, method = "BH"),
    stringsAsFactors = FALSE
  )

  de_results_list[[gse_id]] <- res_df
  write.csv(res_df, sprintf("results/limma_%s.csv", gse_id), row.names = FALSE)
}

saveRDS(de_results_list, "results/limma_5cohorts_list.rds")
capture.output(sessionInfo(), file = "results/sessionInfo_02.txt")
