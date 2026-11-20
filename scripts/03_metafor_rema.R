# ==============================================================================
# 03_metafor_rema.R
# DerSimonian-Laird Random-Effects Meta-Analysis with Knapp-Hartung Adjustment
# Checks finite meta-analysis output for every gene before FDR correction
# ==============================================================================

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
})

de_results_list <- readRDS("results/limma_5cohorts_list.rds")
common_genes <- readRDS("results/common_genes.rds")
cohorts <- names(de_results_list)
stopifnot(length(cohorts) == 5L)

logfc_mat <- sapply(cohorts, function(g) de_results_list[[g]]$logFC)
se_mat <- sapply(cohorts, function(g) de_results_list[[g]]$SE)
rownames(logfc_mat) <- common_genes
rownames(se_mat) <- common_genes

n_genes <- length(common_genes)
meta_rows <- vector("list", n_genes)

for (i in seq_len(n_genes)) {
  yi <- as.numeric(logfc_mat[i, ])
  sei <- as.numeric(se_mat[i, ])
  
  fit <- tryCatch(
    metafor::rma(yi = yi, sei = sei, method = "DL", test = "knha"),
    error = function(e) stop(sprintf("metafor::rma failed on gene %s: %s", common_genes[i], e$message))
  )
  
  # Note: metafor::rma stores the Knapp-Hartung t-statistic (df = 4) in $zval
  meta_rows[[i]] <- data.frame(
    GeneSymbol = common_genes[i],
    k = as.integer(fit$k),
    Pooled_log2FC = as.numeric(fit$beta),
    SE = as.numeric(fit$se),
    tval_knha = as.numeric(fit$zval),
    df = as.integer(fit$ddf),
    pval = as.numeric(fit$pval),
    ci.lb = as.numeric(fit$ci.lb),
    ci.ub = as.numeric(fit$ci.ub),
    tau2 = as.numeric(fit$tau2),
    I2 = as.numeric(fit$I2),
    H2 = as.numeric(fit$H2),
    QE = as.numeric(fit$QE),
    QEp = as.numeric(fit$QEp),
    stringsAsFactors = FALSE
  )
}

all_genes_df <- bind_rows(meta_rows)

# Verify finite estimates prior to Benjamini-Hochberg FDR correction
stopifnot(
  nrow(all_genes_df) == n_genes,
  all(is.finite(all_genes_df$Pooled_log2FC)),
  all(is.finite(all_genes_df$SE) & all_genes_df$SE > 0),
  all(is.finite(all_genes_df$pval) & all_genes_df$pval >= 0 & all_genes_df$pval <= 1),
  all(all_genes_df$df == 4L)
)

all_genes_df$FDR <- p.adjust(all_genes_df$pval, method = "BH")
all_genes_df <- all_genes_df %>% arrange(FDR, desc(abs(Pooled_log2FC)))

consensus_degs_df <- all_genes_df %>%
  filter(FDR < 0.05, abs(Pooled_log2FC) >= 0.50) %>%
  mutate(Direction = ifelse(Pooled_log2FC >= 0.50, "Upregulated", "Downregulated"))

write.csv(all_genes_df, "results/rema_5cohorts_all_genes.csv", row.names = FALSE)
write.csv(consensus_degs_df, "results/rema_5cohorts_consensus_degs.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "results/sessionInfo_03.txt")
