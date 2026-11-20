# ==============================================================================
# 04_clusterProfiler_enrichment.R
# Directional GO:BP Enrichment with Explicit Unique Entrez Universe Reporting
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

all_genes_df <- read.csv("results/rema_5cohorts_all_genes.csv", stringsAsFactors = FALSE)
consensus_degs_df <- read.csv("results/rema_5cohorts_consensus_degs.csv", stringsAsFactors = FALSE)

bg_symbols <- unique(all_genes_df$GeneSymbol)
up_symbols <- unique(consensus_degs_df$GeneSymbol[consensus_degs_df$Pooled_log2FC >= 0.50])
down_symbols <- unique(consensus_degs_df$GeneSymbol[consensus_degs_df$Pooled_log2FC <= -0.50])

map_symbols_to_entrez <- function(symbols) {
  if (length(symbols) == 0L) {
    return(list(
      entrez_ids = character(0),
      mapped_symbols = character(0),
      unmapped_symbols = character(0)
    ))
  }
  mapped_df <- suppressMessages(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  valid_df <- mapped_df[!is.na(mapped_df$SYMBOL) & !is.na(mapped_df$ENTREZID), , drop = FALSE]
  mapped_syms <- unique(as.character(valid_df$SYMBOL))
  unmapped_syms <- setdiff(as.character(symbols), mapped_syms)
  entrez_vec <- unique(as.character(valid_df$ENTREZID))
  list(
    entrez_ids = entrez_vec,
    mapped_symbols = mapped_syms,
    unmapped_symbols = unmapped_syms
  )
}

bg_map <- map_symbols_to_entrez(bg_symbols)
up_map <- map_symbols_to_entrez(up_symbols)
down_map <- map_symbols_to_entrez(down_symbols)

bg_entrez <- bg_map$entrez_ids
up_entrez <- intersect(up_map$entrez_ids, bg_entrez)
down_entrez <- intersect(down_map$entrez_ids, bg_entrez)

n_unmapped_vec <- c(
  length(bg_map$unmapped_symbols),
  length(up_map$unmapped_symbols),
  length(down_map$unmapped_symbols)
)
n_input_vec <- c(length(bg_symbols), length(up_symbols), length(down_symbols))
n_mapped_entrez_vec <- c(length(bg_entrez), length(up_entrez), length(down_entrez))

mapping_audit <- data.frame(
  Gene_Set = c("Eligible_Background_Intersection", "Upregulated_Consensus_DEGs", "Downregulated_Consensus_DEGs"),
  Input_Unique_Symbols = n_input_vec,
  Mapped_Unique_Entrez_Universe = n_mapped_entrez_vec,
  Unmapped_Symbols = n_unmapped_vec,
  Percent_Symbols_Mapped = round(100 * (n_input_vec - n_unmapped_vec) / pmax(n_input_vec, 1L), 2),
  stringsAsFactors = FALSE
)
write.csv(mapping_audit, "results/entrez_mapping_audit.csv", row.names = FALSE)

run_go_bp <- function(entrez_vec, universe_vec, label) {
  if (length(entrez_vec) == 0L || length(universe_vec) == 0L) {
    message(sprintf("[%s] No mapped Entrez IDs available for GO:BP enrichment; writing empty result table.", label))
    return(data.frame(
      ID = character(0), Description = character(0), GeneRatio = character(0),
      BgRatio = character(0), RichFactor = numeric(0), FoldEnrichment = numeric(0),
      oddsRatio = numeric(0), zScore = numeric(0), pvalue = numeric(0),
      p.adjust = numeric(0), qvalue = numeric(0), geneID = character(0),
      Count = integer(0), Status_Note = character(0),
      stringsAsFactors = FALSE
    ))
  }
  ego <- enrichGO(
    gene = entrez_vec,
    universe = universe_vec,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  out_df <- as.data.frame(ego)
  if (nrow(out_df) > 0L) {
    out_df <- out_df %>% filter(p.adjust < 0.05)
  }
  return(out_df)
}

up_go_df <- run_go_bp(up_entrez, bg_entrez, "Upregulated")
down_go_df <- run_go_bp(down_entrez, bg_entrez, "Downregulated")

write.csv(up_go_df, "results/GO_enrichment_upregulated.csv", row.names = FALSE)
write.csv(down_go_df, "results/GO_enrichment_downregulated.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "results/sessionInfo_04.txt")
