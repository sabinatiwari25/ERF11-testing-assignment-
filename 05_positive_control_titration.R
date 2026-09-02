############################################################
## CONTROL 4: POSITIVE CONTROL
## ERF11 Air vs Ethylene RNA-seq
##
## Purpose:
## Plant artificial positive expression changes into genes
## that were originally non-significant, then test whether
## DESeq2 recovers them.
##
## Titration:
## planted log2FC = 0.25, 0.5, 1, 2, 3
##
## Preserve:
## - original genes
## - original samples
## - original biological variation
## - original library structure
## - 3 Air + 3 Ethylene design
## - same filtering
## - same DESeq2 pipeline
##
## Change:
## - counts of selected genes in Ethylene samples only
############################################################


# ==========================================================
# 1. WORKING DIRECTORY
# ==========================================================

setwd(
  "C:/Users/sabin/OneDrive - Auburn University/Computational Biology Coloquium/Assignment 2"
)


# ==========================================================
# 2. PACKAGES
# ==========================================================

library(DESeq2)
library(dplyr)
library(ggplot2)


# ==========================================================
# 3. READ ORIGINAL RAW COUNTS
# ==========================================================

countdata <- as.matrix(
  read.csv(
    "air and ethylene erf11 gene count matrix.csv",
    row.names = "gene_id",
    check.names = FALSE
  )
)


# Clean gene IDs exactly as original analysis
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


# ==========================================================
# 4. READ METADATA
# ==========================================================

coldata <- read.table(
  "Metadata ERF11.txt",
  header = TRUE,
  row.names = 1
)


# Match count matrix to metadata
countdata <- countdata[
  ,
  rownames(coldata),
  drop = FALSE
]


stopifnot(
  identical(
    colnames(countdata),
    rownames(coldata)
  )
)


# Explicitly define reference level
coldata$Treatment <- factor(
  coldata$Treatment,
  levels = c(
    "air_ERF11",
    "ethylene_ERF11"
  )
)


cat("\nSamples per treatment:\n")
print(table(coldata$Treatment))


# ==========================================================
# 5. READ ORIGINAL DESEQ2 RESULTS
# ==========================================================

original_res <- read.csv(
  "1SEPDGESeq_results_ethylene_ERF11_vs_air_ERF11.csv",
  stringsAsFactors = FALSE,
  check.names = TRUE
)


# Rename first column as gene
if ("X" %in% names(original_res)) {
  
  names(original_res)[
    names(original_res) == "X"
  ] <- "gene"
  
} else if (!"gene" %in% names(original_res)) {
  
  names(original_res)[1] <- "gene"
}


# Clean IDs
original_res$gene <- trimws(
  original_res$gene
)

original_res$gene <- gsub(
  "^gene-",
  "",
  original_res$gene
)

original_res$gene <- gsub(
  "\\|.*",
  "",
  original_res$gene
)


# ==========================================================
# 6. SELECT GENES WITH LITTLE ORIGINAL SIGNAL
# ==========================================================
#
# We want genes that:
# - were tested by DESeq2
# - were clearly non-significant
# - originally had very small fold changes
# - have adequate expression
#
# This gives us a relatively neutral baseline into which
# we can plant known effects.
# ==========================================================

candidate_genes <- original_res %>%
  
  filter(
    
    !is.na(padj),
    
    padj > 0.50,
    
    abs(log2FoldChange) < 0.25,
    
    baseMean >= 50,
    
    gene %in% rownames(countdata)
    
  )


# ==========================================================
# 7. REQUIRE ADEQUATE COUNTS IN ALL SIX SAMPLES
# ==========================================================

candidate_min_counts <- apply(
  countdata[
    candidate_genes$gene,
    ,
    drop = FALSE
  ],
  1,
  min
)


candidate_genes$min_count <-
  candidate_min_counts[
    candidate_genes$gene
  ]


candidate_genes <- candidate_genes %>%
  
  filter(
    min_count >= 10
  )


cat(
  "\nSuitable neutral candidate genes:",
  nrow(candidate_genes),
  "\n"
)


if (nrow(candidate_genes) < 20) {
  
  stop(
    "Fewer than 20 suitable genes were found. ",
    "Relax baseMean or min_count criteria slightly."
  )
}


# ==========================================================
# 8. RANDOMLY SELECT 20 GENES
# ==========================================================
#
# Fixed seed ensures reproducibility.
# SAME 20 genes will be used at every effect size.
# ==========================================================

set.seed(20260901)

n_spike <- 20


spike_genes <- sample(
  candidate_genes$gene,
  n_spike,
  replace = FALSE
)


cat("\nSelected planted genes:\n")
print(spike_genes)


# Save baseline characteristics
spike_gene_info <- candidate_genes %>%
  
  filter(
    gene %in% spike_genes
  ) %>%
  
  select(
    gene,
    baseMean,
    log2FoldChange,
    padj,
    min_count
  )


write.csv(
  spike_gene_info,
  "Positive_control_selected_genes.csv",
  row.names = FALSE
)


# ==========================================================
# 9. IDENTIFY ETHYLENE SAMPLES
# ==========================================================

ethylene_samples <- rownames(coldata)[
  coldata$Treatment == "ethylene_ERF11"
]


cat("\nEthylene samples receiving planted signal:\n")
print(ethylene_samples)


stopifnot(
  length(ethylene_samples) == 3
)


# ==========================================================
# 10. DEFINE EFFECT-SIZE TITRATION
# ==========================================================

effect_sizes <- c(
  0.25,
  0.5,
  1,
  2,
  3
)


cat("\nPlanted log2FC values:\n")
print(effect_sizes)


# ==========================================================
# 11. STORAGE FOR RESULTS
# ==========================================================

positive_summary <- data.frame(
  
  planted_log2FC = numeric(),
  
  fold_multiplier = numeric(),
  
  n_planted = integer(),
  
  n_FDR_detected = integer(),
  
  FDR_detection_rate = numeric(),
  
  n_pipeline_detected = integer(),
  
  pipeline_detection_rate = numeric(),
  
  median_estimated_log2FC = numeric(),
  
  median_added_log2FC = numeric()
  
)


positive_gene_results <- list()


# ==========================================================
# 12. RUN POSITIVE-CONTROL TITRATION
# ==========================================================

for (effect in effect_sizes) {
  
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Planting log2FC =",
    effect,
    "\n"
  )
  
  
  # --------------------------------------------------------
  # Copy original count matrix
  # --------------------------------------------------------
  
  spike_counts <- countdata
  
  
  # log2FC -> multiplicative count change
  multiplier <- 2^effect
  
  
  cat(
    "Count multiplier:",
    round(multiplier, 3),
    "\n"
  )
  
  
  # --------------------------------------------------------
  # Plant signal into ETHYLENE samples only
  # --------------------------------------------------------
  
  spike_counts[
    spike_genes,
    ethylene_samples
  ] <- round(
    
    spike_counts[
      spike_genes,
      ethylene_samples
    ] *
      
      multiplier
    
  )
  
  
  # --------------------------------------------------------
  # Create DESeq2 dataset
  # --------------------------------------------------------
  
  dds_spike <- DESeq2::DESeqDataSetFromMatrix(
    
    countData = round(spike_counts),
    
    colData = coldata,
    
    design = ~ Treatment
    
  )
  
  
  # SAME filtering rule as original analysis
  dds_spike <- dds_spike[
    rowSums(
      DESeq2::counts(dds_spike)
    ) > 20,
  ]
  
  
  # --------------------------------------------------------
  # Run SAME DESeq2 pipeline
  # --------------------------------------------------------
  
  dds_spike <- DESeq2::DESeq(
    dds_spike,
    quiet = TRUE
  )
  
  
  res_spike <- DESeq2::results(
    
    dds_spike,
    
    contrast = c(
      "Treatment",
      "ethylene_ERF11",
      "air_ERF11"
    ),
    
    alpha = 0.05
  )
  
  
  res_spike <- as.data.frame(
    res_spike
  )
  
  
  res_spike$gene <- rownames(
    res_spike
  )
  
  
  # --------------------------------------------------------
  # Extract only the planted genes
  # --------------------------------------------------------
  
  planted_result <- res_spike[
    spike_genes,
    ,
    drop = FALSE
  ]
  
  
  planted_result$gene <-
    rownames(planted_result)
  
  
  planted_result$planted_log2FC <-
    effect
  
  
  planted_result$fold_multiplier <-
    multiplier
  
  
  # --------------------------------------------------------
  # Baseline LFC from original analysis
  # --------------------------------------------------------
  
  baseline_lfc <- original_res$log2FoldChange[
    match(
      planted_result$gene,
      original_res$gene
    )
  ]
  
  
  planted_result$baseline_log2FC <-
    baseline_lfc
  
  
  planted_result$added_log2FC <-
    planted_result$log2FoldChange -
    planted_result$baseline_log2FC
  
  
  # --------------------------------------------------------
  # Detection definition 1:
  #
  # Statistically significant in EXPECTED positive direction
  # --------------------------------------------------------
  
  planted_result$FDR_detected <-
    
    !is.na(planted_result$padj) &
    
    planted_result$padj < 0.05 &
    
    planted_result$log2FoldChange > 0
  
  
  # --------------------------------------------------------
  # Detection definition 2:
  #
  # SAME criterion as the real DEG analysis:
  #
  # FDR < 0.05
  # log2FC > 1
  #
  # This is the primary pipeline-recovery measure.
  # --------------------------------------------------------
  
  planted_result$pipeline_detected <-
    
    !is.na(planted_result$padj) &
    
    planted_result$padj < 0.05 &
    
    planted_result$log2FoldChange > 1
  
  
  # --------------------------------------------------------
  # Calculate recovery
  # --------------------------------------------------------
  
  n_FDR <- sum(
    planted_result$FDR_detected,
    na.rm = TRUE
  )
  
  
  n_pipeline <- sum(
    planted_result$pipeline_detected,
    na.rm = TRUE
  )
  
  
  FDR_rate <-
    100 *
    n_FDR /
    n_spike
  
  
  pipeline_rate <-
    100 *
    n_pipeline /
    n_spike
  
  
  median_lfc <- median(
    planted_result$log2FoldChange,
    na.rm = TRUE
  )
  
  
  median_added <- median(
    planted_result$added_log2FC,
    na.rm = TRUE
  )
  
  
  # --------------------------------------------------------
  # Store summary
  # --------------------------------------------------------
  
  positive_summary <- rbind(
    
    positive_summary,
    
    data.frame(
      
      planted_log2FC =
        effect,
      
      fold_multiplier =
        multiplier,
      
      n_planted =
        n_spike,
      
      n_FDR_detected =
        n_FDR,
      
      FDR_detection_rate =
        FDR_rate,
      
      n_pipeline_detected =
        n_pipeline,
      
      pipeline_detection_rate =
        pipeline_rate,
      
      median_estimated_log2FC =
        median_lfc,
      
      median_added_log2FC =
        median_added
      
    )
  )
  
  
  positive_gene_results[[
    as.character(effect)
  ]] <- planted_result
  
  
  cat(
    "FDR-only recovery:",
    n_FDR,
    "/",
    n_spike,
    "=",
    FDR_rate,
    "%\n"
  )
  
  
  cat(
    "Original pipeline recovery:",
    n_pipeline,
    "/",
    n_spike,
    "=",
    pipeline_rate,
    "%\n"
  )
  
  
  cat(
    "Median estimated log2FC:",
    round(median_lfc, 3),
    "\n"
  )
  
  
  cat(
    "Median added log2FC:",
    round(median_added, 3),
    "\n"
  )
  
}


# ==========================================================
# 13. VIEW TITRATION RESULTS
# ==========================================================

cat(
  "\n========================================\n"
)

cat(
  "POSITIVE-CONTROL TITRATION RESULTS\n"
)

cat(
  "========================================\n"
)


print(
  positive_summary
)


# ==========================================================
# 14. COMBINE GENE-LEVEL RESULTS
# ==========================================================

positive_gene_results_df <- bind_rows(
  positive_gene_results
)


# ==========================================================
# 15. SAVE RESULTS
# ==========================================================

write.csv(
  positive_summary,
  "Positive_control_titration_summary.csv",
  row.names = FALSE
)


write.csv(
  positive_gene_results_df,
  "Positive_control_gene_level_results.csv",
  row.names = FALSE
)


# ==========================================================
# 16. POSITIVE-CONTROL DETECTION PLOT
# ==========================================================

p_positive <- ggplot(
  
  positive_summary,
  
  aes(
    x = planted_log2FC,
    y = pipeline_detection_rate
  )
  
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(
      0,
      100,
      20
    )
  ) +
  
  scale_x_continuous(
    breaks = effect_sizes
  ) +
  
  labs(
    
    x = "Planted log2 fold change",
    
    y = "Recovery rate (%)",
    
    title =
      "Positive-control signal titration",
    
    subtitle =
      "Recovery of 20 artificially altered genes using the original DEG criteria",
    
    caption =
      "Detection criterion: FDR < 0.05 and log2FC > 1"
    
  ) +
  
  theme_classic(
    base_size = 13
  )


p_positive


ggsave(
  
  "Positive_control_titration.png",
  
  p_positive,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)


# ==========================================================
# 17. PLOT ESTIMATED VS PLANTED EFFECT
# ==========================================================

p_effect <- ggplot(
  
  positive_summary,
  
  aes(
    x = planted_log2FC,
    y = median_added_log2FC
  )
  
) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  labs(
    
    x = "Planted log2 fold change",
    
    y = "Median recovered added log2FC",
    
    title =
      "Recovery of planted effect sizes",
    
    subtitle =
      "Dashed line indicates perfect recovery"
    
  ) +
  
  theme_classic(
    base_size = 13
  )


p_effect


ggsave(
  
  "Positive_control_effect_recovery.png",
  
  p_effect,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)


# ==========================================================
# 18. SAVE TEXT SUMMARY
# ==========================================================

sink(
  "Positive_control_summary.txt"
)


cat(
  "CONTROL 4: POSITIVE CONTROL\n"
)

cat(
  "===========================\n\n"
)


cat(
  "Purpose:\n"
)

cat(
  "Artificially introduce known expression differences into genes\n"
)

cat(
  "that originally showed little evidence of an Air/Ethylene effect.\n\n"
)


cat(
  "Signal planted:\n"
)

cat(
  "Counts for the same 20 genes were multiplied in all three\n"
)

cat(
  "ethylene samples by factors corresponding to planted log2FC\n"
)

cat(
  "values of 0.25, 0.5, 1, 2, and 3.\n\n"
)


cat(
  "Preserved:\n"
)

cat(
  "- Original raw count matrix structure\n"
)

cat(
  "- Six samples\n"
)

cat(
  "- Three Air and three Ethylene samples\n"
)

cat(
  "- Original biological variability\n"
)

cat(
  "- Gene identities\n"
)

cat(
  "- DESeq2 analysis method\n"
)

cat(
  "- Filtering rule\n"
)

cat(
  "- Original DEG criteria\n\n"
)


cat(
  "Selected planted genes:",
  n_spike,
  "\n\n"
)


cat(
  "Titration results:\n"
)

print(
  positive_summary
)


cat(
  "\nPrimary recovery criterion:\n"
)

cat(
  "FDR < 0.05 and log2FoldChange > 1\n"
)


sink()


# ==========================================================
# 19. FINAL CONSOLE OUTPUT
# ==========================================================

cat(
  "\n========================================\n"
)

cat(
  "POSITIVE CONTROL COMPLETE\n"
)

cat(
  "========================================\n"
)


print(
  positive_summary
)


cat(
  "\nFiles created:\n"
)

cat(
  "1. Positive_control_selected_genes.csv\n"
)

cat(
  "2. Positive_control_titration_summary.csv\n"
)

cat(
  "3. Positive_control_gene_level_results.csv\n"
)

cat(
  "4. Positive_control_titration.png\n"
)

cat(
  "5. Positive_control_effect_recovery.png\n"
)

cat(
  "6. Positive_control_summary.txt\n"
)

positive_summary
