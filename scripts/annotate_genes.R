# ================================================
# Gene Annotation Script
# Tool: biomaRt (Ensembl)
# Input: results/de_analysis/de_results.csv
# Output: results/de_analysis/de_results_annotated.csv
# ================================================

library(biomaRt)

cat("=================================\n")
cat("      GENE ANNOTATION\n")
cat("=================================\n")

# ── Load DE results ───────────────────────────
de_results <- read.csv("results/de_analysis/de_results.csv", 
                        row.names=1)
colnames(de_results) <- c("baseMean", "log2FoldChange", 
                           "lfcSE", "stat", "pvalue", "padj")

cat(sprintf("Loaded %d transcripts\n", nrow(de_results)))

# ── Clean transcript IDs (remove version numbers) ─
transcript_ids <- rownames(de_results)
transcript_ids_clean <- gsub("\\..*", "", transcript_ids)

# ── Connect to Ensembl ────────────────────────
cat("Connecting to Ensembl biomaRt...\n")
mart <- useEnsembl(biomart="ensembl", 
                   dataset="hsapiens_gene_ensembl",
                   mirror="useast")

# ── Query gene annotations ────────────────────
cat("Querying gene annotations...\n")
annotations <- getBM(
  attributes=c("ensembl_transcript_id",
               "hgnc_symbol",
               "description",
               "gene_biotype",
               "chromosome_name"),
  filters="ensembl_transcript_id",
  values=transcript_ids_clean,
  mart=mart
)

cat(sprintf("Annotations retrieved for %d transcripts\n", 
            nrow(annotations)))

# ── Merge with DE results ─────────────────────
de_results$transcript_id_clean <- transcript_ids_clean

merged <- merge(de_results, annotations,
                by.x="transcript_id_clean",
                by.y="ensembl_transcript_id",
                all.x=TRUE)

rownames(merged) <- transcript_ids[match(merged$transcript_id_clean, 
                                          transcript_ids_clean)]

# ── Save annotated results ────────────────────
write.csv(merged, "results/de_analysis/de_results_annotated.csv",
          row.names=TRUE)

cat("Saved: results/de_analysis/de_results_annotated.csv\n")

# ── Print top 10 most significant with gene names ─
cat("\nTop 10 most significant DE genes:\n")
cat("=================================\n")
sig <- merged[!is.na(merged$padj), ]
sig <- sig[order(sig$padj), ]
top10 <- head(sig[, c("hgnc_symbol", "log2FoldChange", "padj", "description")], 10)
print(top10)

cat("\nDone! ✅\n")

# ── Save summary files ────────────────────────
cat("\nSaving summary files...\n")

# Filter significant genes
sig <- merged[!is.na(merged$padj) & merged$padj < 0.05, ]
sig <- sig[order(sig$padj), ]

# Top 20 overall
top20 <- head(sig[, c("hgnc_symbol", "log2FoldChange", "padj", "description")], 20)
write.csv(top20, "results/de_analysis/top20_significant.csv", row.names=TRUE)
cat("Saved: top20_significant.csv\n")

# Top upregulated (LFC > 1)
up <- sig[sig$log2FoldChange > 1, c("hgnc_symbol", "log2FoldChange", "padj", "description")]
write.csv(head(up, 20), "results/de_analysis/top_upregulated.csv", row.names=TRUE)
cat(sprintf("Saved: top_upregulated.csv (%d genes)\n", nrow(up)))

# Top downregulated (LFC < -1)
dn <- sig[sig$log2FoldChange < -1, c("hgnc_symbol", "log2FoldChange", "padj", "description")]
write.csv(head(dn, 20), "results/de_analysis/top_downregulated.csv", row.names=TRUE)
cat(sprintf("Saved: top_downregulated.csv (%d genes)\n", nrow(dn)))

cat("\nAll files saved! ✅\n")
