############################################################
## CONTROL 5: REDUNDANCY
## FINAL edgeR vs DESeq2 ANALYSIS
##
## Biological comparison:
## Ethylene ERF11 vs Air ERF11
##
## WHY ARE WE DOING THIS?
##
## The primary differential-expression analysis was performed
## with DESeq2. For a redundancy control, we analyze the SAME
## biological experiment using a second established RNA-seq
## method: edgeR.
##
## The purpose is NOT to force edgeR to reproduce exactly the
## same list of significant genes as DESeq2.
##
## Instead, we ask:
##
## 1. Do the two methods estimate similar log2 fold changes?
## 2. Do they agree on the direction of expression changes?
## 3. How much do their statistically significant DEG lists
##    overlap?
##
## IMPORTANT:
##
## Here edgeR is allowed to use its recommended filterByExpr()
## filtering procedure. This is preferable for the final edgeR
## analysis because lowly expressed genes can contribute little
## statistical information while increasing the multiple-testing
## burden.
##
## Shared between the two approaches:
## - same raw count matrix
## - same six biological samples
## - same Air/Ethylene labels
## - same biological comparison
##
## Different:
## - normalization
## - expression filtering
## - dispersion estimation
## - statistical testing
##
## Therefore, this is a REDUNDANCY / SENSITIVITY ANALYSIS,
## not a perfectly independent replication.
############################################################


# ==========================================================
# 1. WORKING DIRECTORY
# ==========================================================

setwd(
  "C:/Users/sabin/OneDrive - Auburn University/Computational Biology Coloquium/Assignment 2"
)


# ==========================================================
# 2. LOAD PACKAGES
# ==========================================================
#
# edgeR:
#   performs the second differential-expression analysis.
#
# dplyr:
#   used for organizing and comparing results.
#
# ggplot2:
#   used for redundancy figures.
# ==========================================================

# Run ONLY if edgeR is not already installed:
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
BiocManager::install("edgeR")

library(edgeR)
library(dplyr)
library(ggplot2)


cat("\n========================================\n")
cat("CONTROL 5: REDUNDANCY - edgeR vs DESeq2\n")
cat("========================================\n")


# ==========================================================
# 3. READ ORIGINAL RAW COUNT MATRIX
# ==========================================================
#
# WHY?
#
# edgeR, like DESeq2, should start with RAW integer counts.
# We do NOT use RPKM or other previously normalized data.
# ==========================================================

countdata <- as.matrix(
  read.csv(
    "air and ethylene erf11 gene count matrix.csv",
    row.names = "gene_id",
    check.names = FALSE
  )
)


# Clean gene IDs exactly as in the DESeq2 analysis
rownames(countdata) <- gsub(
  "^gene-",
  "",
  rownames(countdata)
)

rownames(countdata) <- gsub(
  "\\|.*",
  "",
  rownames(countdata)
)


storage.mode(countdata) <- "numeric"


cat(
  "\nGenes in original raw-count matrix:",
  nrow(countdata),
  "\n"
)


# ==========================================================
# 4. READ METADATA
# ==========================================================

coldata <- read.table(
  "Metadata ERF11.txt",
  header = TRUE,
  row.names = 1
)


# Match count columns to metadata rows
countdata <- countdata[
  ,
  rownames(coldata),
  drop = FALSE
]


# ==========================================================
# 5. VERIFY SAMPLE MATCHING
# ==========================================================
#
# WHY?
#
# A differential-expression analysis is meaningless if the
# count columns and metadata rows are assigned to different
# biological samples.
# ==========================================================

stopifnot(
  identical(
    colnames(countdata),
    rownames(coldata)
  )
)


cat(
  "Count matrix and metadata order identical:",
  identical(
    colnames(countdata),
    rownames(coldata)
  ),
  "\n"
)


# ==========================================================
# 6. DEFINE TREATMENT REFERENCE
# ==========================================================
#
# Air is the reference.
#
# Therefore:
#
# POSITIVE logFC = higher expression in ETHYLENE
# NEGATIVE logFC = lower expression in ETHYLENE
# ==========================================================

coldata$Treatment <- factor(
  coldata$Treatment,
  levels = c(
    "air_ERF11",
    "ethylene_ERF11"
  )
)


cat("\nSamples per treatment:\n")

print(
  table(coldata$Treatment)
)


stopifnot(
  sum(coldata$Treatment == "air_ERF11") == 3,
  sum(coldata$Treatment == "ethylene_ERF11") == 3
)


# ==========================================================
# 7. CREATE edgeR DGEList
# ==========================================================
#
# DGEList stores:
#
# - raw gene counts
# - sample information
# - library sizes
# - normalization information
# ==========================================================

y <- edgeR::DGEList(
  counts = round(countdata),
  group = coldata$Treatment
)


cat("\nOriginal library sizes:\n")

print(
  y$samples[, c("group", "lib.size")]
)


# ==========================================================
# 8. CREATE DESIGN MATRIX
# ==========================================================
#
# The design tells edgeR that expression is modeled as a
# function of Treatment.
# ==========================================================

design <- model.matrix(
  ~ Treatment,
  data = coldata
)


cat("\nDesign matrix:\n")

print(design)


cat("\nDesign coefficients:\n")

print(
  colnames(design)
)


# Check that the expected coefficient exists
stopifnot(
  "Treatmentethylene_ERF11" %in%
    colnames(design)
)


# ==========================================================
# 9. edgeR-SPECIFIC EXPRESSION FILTER
# ==========================================================
#
# WHY filterByExpr()?
#
# Genes with extremely low expression contain little power
# for differential-expression testing.
#
# Testing thousands of such genes:
#
# - increases the multiple-testing burden
# - may reduce statistical power
# - can complicate dispersion estimation
#
# edgeR provides filterByExpr() specifically to determine
# whether expression is sufficiently high for the experimental
# design.
#
# IMPORTANT:
#
# We are NOT doing this to make edgeR agree with DESeq2.
# We are using edgeR according to its recommended workflow.
# ==========================================================

keep_edgeR <- edgeR::filterByExpr(
  y,
  design = design
)


cat(
  "\nGenes before edgeR filtering:",
  nrow(y),
  "\n"
)


cat(
  "Genes retained by filterByExpr():",
  sum(keep_edgeR),
  "\n"
)


cat(
  "Genes removed by filterByExpr():",
  sum(!keep_edgeR),
  "\n"
)


# Keep sufficiently expressed genes
y <- y[
  keep_edgeR,
  ,
  keep.lib.sizes = FALSE
]


# ==========================================================
# 10. TMM NORMALIZATION
# ==========================================================
#
# WHY?
#
# Sequencing libraries have different depths.
#
# TMM = Trimmed Mean of M-values.
#
# It adjusts effective library sizes so that expression
# differences are not simply caused by differences in
# sequencing depth or composition.
# ==========================================================

y <- edgeR::calcNormFactors(
  y,
  method = "TMM"
)


cat("\nedgeR sample normalization information:\n")

print(
  y$samples
)


# ==========================================================
# 11. ESTIMATE DISPERSION
# ==========================================================
#
# WHY?
#
# Biological replicates are not identical.
#
# edgeR models counts with a negative-binomial distribution.
# Dispersion quantifies biological + technical variation above
# simple Poisson counting variation.
#
# With only 3 biological replicates per treatment, dispersion
# estimation is especially important because large variability
# reduces statistical confidence.
# ==========================================================

y <- edgeR::estimateDisp(
  y,
  design
)


cat(
  "\nCommon dispersion:",
  y$common.dispersion,
  "\n"
)


common_BCV <- sqrt(
  y$common.dispersion
)


cat(
  "Common biological coefficient of variation (BCV):",
  common_BCV,
  "\n"
)


# Save BCV plot
png(
  "edgeR_BCV_plot.png",
  width = 1800,
  height = 1500,
  res = 250
)

plotBCV(y)

dev.off()


# ==========================================================
# 12. FIT QUASI-LIKELIHOOD MODEL
# ==========================================================
#
# WHY QL?
#
# edgeR's quasi-likelihood framework accounts for uncertainty
# in gene-specific variability.
#
# glmQLFit() + glmQLFTest() is generally preferred for
# differential-expression inference because it provides
# rigorous control of false positives.
#
# robust = TRUE provides additional protection against genes
# with unusually large dispersions.
# ==========================================================

fit <- edgeR::glmQLFit(
  y,
  design,
  robust = TRUE
)


# ==========================================================
# 13. TEST ETHYLENE VS AIR
# ==========================================================
#
# Because Air is the reference level:
#
# Treatmentethylene_ERF11 > 0
#     means higher expression under ethylene.
#
# Treatmentethylene_ERF11 < 0
#     means lower expression under ethylene.
# ==========================================================

qlf <- edgeR::glmQLFTest(
  fit,
  coef = "Treatmentethylene_ERF11"
)


# ==========================================================
# 14. EXTRACT ALL edgeR RESULTS
# ==========================================================

edgeR_results <- edgeR::topTags(
  qlf,
  n = Inf,
  sort.by = "none"
)$table


edgeR_results$gene <- rownames(
  edgeR_results
)


# Rename using base R to avoid package masking problems
names(edgeR_results)[
  names(edgeR_results) == "logFC"
] <- "edgeR_log2FC"


names(edgeR_results)[
  names(edgeR_results) == "PValue"
] <- "edgeR_PValue"


names(edgeR_results)[
  names(edgeR_results) == "FDR"
] <- "edgeR_FDR"


cat("\nedgeR result columns:\n")

print(
  names(edgeR_results)
)


# ==========================================================
# 15. DEFINE SIGNIFICANCE CRITERIA
# ==========================================================
#
# We use the SAME final biological criterion used in the
# DESeq2 analysis:
#
# FDR < 0.05
#
# AND
#
# |log2FC| > 1
# ==========================================================

alpha <- 0.05
lfc_cut <- 1


edgeR_results <- edgeR_results %>%
  
  dplyr::mutate(
    
    edgeR_group = dplyr::case_when(
      
      !is.na(edgeR_FDR) &
        edgeR_FDR < alpha &
        edgeR_log2FC > lfc_cut ~
        "up",
      
      !is.na(edgeR_FDR) &
        edgeR_FDR < alpha &
        edgeR_log2FC < -lfc_cut ~
        "down",
      
      TRUE ~
        "ns"
    )
  )


# ==========================================================
# 16. COUNT edgeR DEGs
# ==========================================================

n_edgeR_up <- sum(
  edgeR_results$edgeR_group == "up",
  na.rm = TRUE
)


n_edgeR_down <- sum(
  edgeR_results$edgeR_group == "down",
  na.rm = TRUE
)


n_edgeR_DEG <-
  n_edgeR_up +
  n_edgeR_down


cat("\n========================================\n")
cat("STANDARD edgeR RESULTS\n")
cat("========================================\n")


cat(
  "Genes tested:",
  nrow(edgeR_results),
  "\n"
)


cat(
  "Upregulated:",
  n_edgeR_up,
  "\n"
)


cat(
  "Downregulated:",
  n_edgeR_down,
  "\n"
)


cat(
  "Total edgeR DEGs:",
  n_edgeR_DEG,
  "\n"
)


# ==========================================================
# 17. ALSO COUNT FDR-ONLY RESULTS
# ==========================================================
#
# WHY?
#
# This tells us whether genes fail because of:
#
# - statistical significance itself
#
# or
#
# - the additional |log2FC| > 1 threshold.
# ==========================================================

n_edgeR_FDR_only <- sum(
  edgeR_results$edgeR_FDR < 0.05,
  na.rm = TRUE
)


cat(
  "Genes with edgeR FDR < 0.05 regardless of fold change:",
  n_edgeR_FDR_only,
  "\n"
)


# ==========================================================
# 18. LOOK AT STRONGEST edgeR RESULTS
# ==========================================================
#
# WHY?
#
# If edgeR finds few FDR-significant genes, inspecting the
# strongest results tells us whether many genes still have
# low nominal P-values but fail after multiple-testing
# correction.
# ==========================================================

top_edgeR <- edgeR::topTags(
  qlf,
  n = 20
)$table


cat("\nTop 20 edgeR results:\n")

print(
  top_edgeR
)


write.csv(
  top_edgeR,
  "edgeR_top20_results.csv",
  row.names = TRUE
)


# ==========================================================
# 19. SAVE FULL edgeR RESULTS
# ==========================================================

write.csv(
  edgeR_results,
  "edgeR_filterByExpr_ERF11_ethylene_vs_air_results.csv",
  row.names = FALSE
)


# ==========================================================
# 20. READ ORIGINAL DESeq2 RESULTS
# ==========================================================
#
# These are the results from the primary analysis.
# ==========================================================

deseq <- read.csv(
  "1SEPDGESeq_results_ethylene_ERF11_vs_air_ERF11.csv",
  stringsAsFactors = FALSE,
  check.names = TRUE
)


# Identify gene-ID column
if ("X" %in% names(deseq)) {
  
  names(deseq)[
    names(deseq) == "X"
  ] <- "gene"
  
} else if (!"gene" %in% names(deseq)) {
  
  names(deseq)[1] <- "gene"
}


# Same gene-ID cleaning
deseq$gene <- trimws(
  deseq$gene
)


deseq$gene <- gsub(
  "^gene-",
  "",
  deseq$gene
)


deseq$gene <- gsub(
  "\\|.*",
  "",
  deseq$gene
)


# Rename columns safely
names(deseq)[
  names(deseq) == "log2FoldChange"
] <- "DESeq2_log2FC"


names(deseq)[
  names(deseq) == "padj"
] <- "DESeq2_padj"


names(deseq)[
  names(deseq) == "pvalue"
] <- "DESeq2_pvalue"


# ==========================================================
# 21. DEFINE DESeq2 DEG CLASSIFICATION
# ==========================================================

deseq <- deseq %>%
  
  dplyr::mutate(
    
    DESeq2_group = dplyr::case_when(
      
      !is.na(DESeq2_padj) &
        DESeq2_padj < alpha &
        DESeq2_log2FC > lfc_cut ~
        "up",
      
      !is.na(DESeq2_padj) &
        DESeq2_padj < alpha &
        DESeq2_log2FC < -lfc_cut ~
        "down",
      
      TRUE ~
        "ns"
    )
  )


# ==========================================================
# 22. MERGE ONLY GENES TESTED BY BOTH METHODS
# ==========================================================
#
# WHY?
#
# edgeR filterByExpr() removes low-information genes.
#
# Therefore comparisons should be made only among genes that
# were actually available in both analyses.
# ==========================================================

comparison <- dplyr::inner_join(
  
  deseq %>%
    
    dplyr::select(
      gene,
      DESeq2_log2FC,
      DESeq2_padj,
      DESeq2_group
    ),
  
  edgeR_results %>%
    
    dplyr::select(
      gene,
      edgeR_log2FC,
      edgeR_FDR,
      edgeR_group
    ),
  
  by = "gene"
)


cat(
  "\nGenes tested by both DESeq2 and edgeR:",
  nrow(comparison),
  "\n"
)


# ==========================================================
# 23. REMOVE NON-FINITE VALUES FOR CORRELATION
# ==========================================================

comparison_complete <- comparison %>%
  
  dplyr::filter(
    
    is.finite(DESeq2_log2FC),
    
    is.finite(edgeR_log2FC)
  )


# ==========================================================
# 24. CORRELATION ACROSS ALL COMMON GENES
# ==========================================================
#
# WHY?
#
# DEG cutoffs can be sensitive to FDR and dispersion.
#
# Fold-change correlation tests a more fundamental question:
#
# Do both methods estimate approximately the same biological
# expression response?
# ==========================================================

pearson_cor <- cor(
  
  comparison_complete$DESeq2_log2FC,
  
  comparison_complete$edgeR_log2FC,
  
  method = "pearson"
)


spearman_cor <- cor(
  
  comparison_complete$DESeq2_log2FC,
  
  comparison_complete$edgeR_log2FC,
  
  method = "spearman"
)


cat(
  "\nPearson log2FC correlation:",
  round(pearson_cor, 4),
  "\n"
)


cat(
  "Spearman log2FC correlation:",
  round(spearman_cor, 4),
  "\n"
)


# ==========================================================
# 25. OVERALL DIRECTION AGREEMENT
# ==========================================================
#
# Example:
#
# DESeq2 = positive
# edgeR  = positive
#
# -> same biological direction
# ==========================================================

overall_direction_agreement <- mean(
  
  sign(
    comparison_complete$DESeq2_log2FC
  ) ==
    
    sign(
      comparison_complete$edgeR_log2FC
    )
  
) * 100


cat(
  "Overall log2FC direction agreement:",
  round(overall_direction_agreement, 1),
  "%\n"
)


# ==========================================================
# 26. DESeq2 SIGNIFICANT GENES AMONG COMMON GENE SET
# ==========================================================

deseq_sig_compare <- comparison_complete %>%
  
  dplyr::filter(
    DESeq2_group != "ns"
  )


n_deseq_common_DEGs <- nrow(
  deseq_sig_compare
)


cat(
  "\nDESeq2 DEGs retained/tested by edgeR:",
  n_deseq_common_DEGs,
  "\n"
)


# ==========================================================
# 27. DIRECTION AGREEMENT FOR DESeq2 DEGs
# ==========================================================
#
# This is especially important.
#
# edgeR may not call a DESeq2 DEG statistically significant,
# but it may still estimate the SAME biological direction.
# ==========================================================

direction_DESeq2_DEGs <- mean(
  
  sign(
    deseq_sig_compare$DESeq2_log2FC
  ) ==
    
    sign(
      deseq_sig_compare$edgeR_log2FC
    )
  
) * 100


cat(
  "Direction agreement among DESeq2 DEGs:",
  round(direction_DESeq2_DEGs, 1),
  "%\n"
)


# ==========================================================
# 28. CORRELATION AMONG DESeq2 DEGs
# ==========================================================

pearson_DESeq2_DEGs <- cor(
  
  deseq_sig_compare$DESeq2_log2FC,
  
  deseq_sig_compare$edgeR_log2FC,
  
  method = "pearson"
)


spearman_DESeq2_DEGs <- cor(
  
  deseq_sig_compare$DESeq2_log2FC,
  
  deseq_sig_compare$edgeR_log2FC,
  
  method = "spearman"
)


cat(
  "Pearson correlation among DESeq2 DEGs:",
  round(pearson_DESeq2_DEGs, 4),
  "\n"
)


cat(
  "Spearman correlation among DESeq2 DEGs:",
  round(spearman_DESeq2_DEGs, 4),
  "\n"
)


# ==========================================================
# 29. COUNT SIGNIFICANT GENES IN BOTH METHODS
# ==========================================================

DESeq2_gene_set <- comparison$gene[
  comparison$DESeq2_group != "ns"
]


edgeR_gene_set <- comparison$gene[
  comparison$edgeR_group != "ns"
]


shared_genes <- intersect(
  DESeq2_gene_set,
  edgeR_gene_set
)


n_DESeq2 <- length(
  DESeq2_gene_set
)


n_edgeR <- length(
  edgeR_gene_set
)


n_shared <- length(
  shared_genes
)


cat("\n========================================\n")
cat("DEG OVERLAP\n")
cat("========================================\n")


cat(
  "DESeq2 DEGs:",
  n_DESeq2,
  "\n"
)


cat(
  "edgeR DEGs:",
  n_edgeR,
  "\n"
)


cat(
  "Shared DEGs:",
  n_shared,
  "\n"
)


# ==========================================================
# 30. JACCARD OVERLAP
# ==========================================================

union_DEGs <- union(
  DESeq2_gene_set,
  edgeR_gene_set
)


if (length(union_DEGs) > 0) {
  
  jaccard <- n_shared /
    length(union_DEGs)
  
} else {
  
  jaccard <- NA_real_
}


cat(
  "Jaccard overlap:",
  round(jaccard, 4),
  "\n"
)


# ==========================================================
# 31. SIGNIFICANCE CLASSIFICATION TABLE
# ==========================================================

classification_table <- table(
  
  DESeq2 =
    comparison$DESeq2_group,
  
  edgeR =
    comparison$edgeR_group
)


cat(
  "\nClassification table:\n"
)


print(
  classification_table
)


# ==========================================================
# 32. SAVE COMPARISON TABLE
# ==========================================================

write.csv(
  comparison,
  "DESeq2_vs_edgeR_filterByExpr_comparison.csv",
  row.names = FALSE
)


# ==========================================================
# 33. PLOT DESeq2 vs edgeR LOG2FC
# ==========================================================

p_logfc <- ggplot(
  
  comparison_complete,
  
  aes(
    x = DESeq2_log2FC,
    y = edgeR_log2FC
  )
  
) +
  
  geom_point(
    alpha = 0.3,
    size = 1.2
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    
    x =
      "DESeq2 log2 fold change",
    
    y =
      "edgeR log2 fold change",
    
    title =
      "Redundancy control: DESeq2 vs edgeR",
    
    subtitle =
      paste0(
        "Pearson r = ",
        round(pearson_cor, 3),
        "; Spearman rho = ",
        round(spearman_cor, 3)
      ),
    
    caption =
      "edgeR uses filterByExpr, TMM normalization and quasi-likelihood testing."
    
  ) +
  
  theme_classic(
    base_size = 13
  )


p_logfc


ggsave(
  "Redundancy_FINAL_DESeq2_vs_edgeR_log2FC.png",
  p_logfc,
  width = 8,
  height = 7,
  dpi = 300
)


# ==========================================================
# 34. PLOT ONLY THE DESeq2 SIGNIFICANT GENES
# ==========================================================
#
# This shows whether edgeR sees the same direction/magnitude
# for genes called significant by DESeq2.
# ==========================================================

p_deseq_DEGs <- ggplot(
  
  deseq_sig_compare,
  
  aes(
    x = DESeq2_log2FC,
    y = edgeR_log2FC
  )
  
) +
  
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    
    x =
      "DESeq2 log2 fold change",
    
    y =
      "edgeR log2 fold change",
    
    title =
      "edgeR effect estimates for DESeq2 significant genes",
    
    subtitle =
      paste0(
        "Direction agreement = ",
        round(direction_DESeq2_DEGs, 1),
        "%"
      ),
    
    caption =
      "Genes shown passed the DESeq2 DEG criterion and were retained by edgeR."
    
  ) +
  
  theme_classic(
    base_size = 13
  )


p_deseq_DEGs


ggsave(
  "Redundancy_FINAL_DESeq2_DEGs_edgeR_effects.png",
  p_deseq_DEGs,
  width = 8,
  height = 7,
  dpi = 300
)


# ==========================================================
# 35. SAVE FINAL SUMMARY
# ==========================================================

sink(
  "Redundancy_FINAL_summary.txt"
)


cat(
  "CONTROL 5: REDUNDANCY\n"
)

cat(
  "=====================\n\n"
)


cat(
  "Biological comparison:\n"
)

cat(
  "Ethylene ERF11 vs Air ERF11\n\n"
)


cat(
  "Purpose:\n"
)

cat(
  "Assess whether differential-expression effect estimates are robust\n"
)

cat(
  "to analysis with a second established RNA-seq method.\n\n"
)


cat(
  "Primary method:\n"
)

cat(
  "DESeq2\n\n"
)


cat(
  "Redundant method:\n"
)

cat(
  "edgeR using filterByExpr, TMM normalization, dispersion estimation,\n"
)

cat(
  "and quasi-likelihood testing.\n\n"
)


cat(
  "Why filterByExpr was used:\n"
)

cat(
  "edgeR recommends filtering genes with insufficient expression before\n"
)

cat(
  "dispersion estimation and multiple-testing correction. This reduces\n"
)

cat(
  "the contribution of genes with little statistical information.\n\n"
)


cat(
  "edgeR filtering:\n"
)

cat(
  "Original genes:",
  nrow(countdata),
  "\n"
)

cat(
  "Genes retained:",
  sum(keep_edgeR),
  "\n"
)

cat(
  "Genes removed:",
  sum(!keep_edgeR),
  "\n\n"
)


cat(
  "edgeR dispersion:\n"
)

cat(
  "Common dispersion:",
  y$common.dispersion,
  "\n"
)

cat(
  "Common BCV:",
  common_BCV,
  "\n\n"
)


cat(
  "Differential-expression results:\n"
)

cat(
  "DESeq2 DEGs among common genes:",
  n_DESeq2,
  "\n"
)

cat(
  "edgeR DEGs:",
  n_edgeR,
  "\n"
)

cat(
  "Shared significant DEGs:",
  n_shared,
  "\n"
)

cat(
  "Jaccard overlap:",
  round(jaccard, 4),
  "\n\n"
)


cat(
  "Fold-change agreement across common genes:\n"
)

cat(
  "Pearson correlation:",
  round(pearson_cor, 4),
  "\n"
)

cat(
  "Spearman correlation:",
  round(spearman_cor, 4),
  "\n"
)

cat(
  "Overall directional agreement:",
  round(overall_direction_agreement, 1),
  "%\n\n"
)


cat(
  "Agreement among DESeq2 significant genes:\n"
)

cat(
  "Direction agreement:",
  round(direction_DESeq2_DEGs, 1),
  "%\n"
)

cat(
  "Pearson correlation:",
  round(pearson_DESeq2_DEGs, 4),
  "\n"
)

cat(
  "Spearman correlation:",
  round(spearman_DESeq2_DEGs, 4),
  "\n\n"
)


cat(
  "Important limitation:\n"
)

cat(
  "DESeq2 and edgeR are not fully independent because they share the\n"
)

cat(
  "same raw counts, biological samples, treatment labels and general\n"
)

cat(
  "negative-binomial framework. Upstream errors such as incorrect sample\n"
)

cat(
  "labels could therefore affect both approaches in the same way.\n\n"
)


cat(
  "Classification table:\n"
)

print(
  classification_table
)


sink()


# ==========================================================
# 36. FINAL CONSOLE OUTPUT
# ==========================================================

cat("\n========================================\n")
cat("FINAL REDUNDANCY ANALYSIS COMPLETE\n")
cat("========================================\n")


cat(
  "Genes retained by filterByExpr:",
  sum(keep_edgeR),
  "\n"
)


cat(
  "Common BCV:",
  round(common_BCV, 4),
  "\n"
)


cat(
  "DESeq2 DEGs:",
  n_DESeq2,
  "\n"
)


cat(
  "edgeR DEGs:",
  n_edgeR,
  "\n"
)


cat(
  "Shared DEGs:",
  n_shared,
  "\n"
)


cat(
  "Pearson correlation:",
  round(pearson_cor, 4),
  "\n"
)


cat(
  "Spearman correlation:",
  round(spearman_cor, 4),
  "\n"
)


cat(
  "Overall direction agreement:",
  round(overall_direction_agreement, 1),
  "%\n"
)


cat(
  "Direction agreement among DESeq2 DEGs:",
  round(direction_DESeq2_DEGs, 1),
  "%\n"
)


cat(
  "Pearson correlation among DESeq2 DEGs:",
  round(pearson_DESeq2_DEGs, 4),
  "\n"
)


cat(
  "Spearman correlation among DESeq2 DEGs:",
  round(spearman_DESeq2_DEGs, 4),
  "\n"
)


cat("\nFiles created:\n")

cat(
  "1. edgeR_filterByExpr_ERF11_ethylene_vs_air_results.csv\n"
)

cat(
  "2. edgeR_top20_results.csv\n"
)

cat(
  "3. edgeR_BCV_plot.png\n"
)

cat(
  "4. DESeq2_vs_edgeR_filterByExpr_comparison.csv\n"
)

cat(
  "5. Redundancy_FINAL_DESeq2_vs_edgeR_log2FC.png\n"
)

cat(
  "6. Redundancy_FINAL_DESeq2_DEGs_edgeR_effects.png\n"
)

cat(
  "7. Redundancy_FINAL_summary.txt\n"
)

