############################################################
## WEEK 2 ASSIGNMENT: TESTING MINDSET
##
## MAIN ANALYSIS + CONTROL 1: KNOWN ANSWER
##
## Comparison:
## Ethylene ERF11 vs Air ERF11
##
## Known-answer control:
## Independently identified EIN3-regulated genes
##
## DEG definition:
## FDR < 0.05
## |log2FoldChange| > 1
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

library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)


# ==========================================================
# 3. READ DESEQ2 RESULTS
# ==========================================================

df <- read.csv(
  "1SEPDGESeq_results_ethylene_ERF11_vs_air_ERF11.csv",
  stringsAsFactors = FALSE,
  check.names = TRUE
)


# Check imported column names
cat("\nDESeq2 result columns:\n")
print(names(df))


# ----------------------------------------------------------
# Identify gene-ID column
# ----------------------------------------------------------

# If row names were exported by write.csv(), the first
# column is usually called X.

if ("X" %in% names(df)) {
  
  names(df)[names(df) == "X"] <- "gene"
  
} else if (!"gene" %in% names(df)) {
  
  names(df)[1] <- "gene"
}


# Clean gene IDs in DESeq2 results
df$gene <- trimws(df$gene)

df$gene <- gsub(
  "^gene-",
  "",
  df$gene
)

df$gene <- gsub(
  "\\|.*",
  "",
  df$gene
)


cat("\nFirst DESeq2 gene IDs:\n")
print(head(df$gene))


# ==========================================================
# 4. READ EIN3-REGULATED GENE LIST
# ==========================================================

genes_to_label <- read.csv(
  "EIN3_ethylene.csv",
  stringsAsFactors = FALSE,
  check.names = TRUE
)


cat("\nEIN3 file columns:\n")
print(names(genes_to_label))


# ----------------------------------------------------------
# IMPORTANT:
# Keep ONLY the two columns needed for the analysis.
#
# This fixes the error:
# "Input columns in y must be unique"
# ----------------------------------------------------------

required_cols <- c(
  "Gene_ID",
  "Gene_Name"
)


if (!all(required_cols %in% names(genes_to_label))) {
  
  stop(
    "EIN3_ethylene.csv must contain Gene_ID and Gene_Name columns."
  )
}


genes_to_label <- genes_to_label %>%
  
  select(
    Gene_ID,
    Gene_Name
  )


# Clean EIN3 annotations
genes_to_label <- genes_to_label %>%
  
  mutate(
    
    Gene_ID = trimws(Gene_ID),
    
    Gene_Name = trimws(Gene_Name)
    
  )


# Apply same gene-ID cleanup as count matrix / DESeq2
genes_to_label$Gene_ID <- gsub(
  "^gene-",
  "",
  genes_to_label$Gene_ID
)

genes_to_label$Gene_ID <- gsub(
  "\\|.*",
  "",
  genes_to_label$Gene_ID
)


# Remove missing IDs
genes_to_label <- genes_to_label %>%
  
  filter(
    
    !is.na(Gene_ID),
    
    Gene_ID != ""
    
  )


# Keep only one annotation per Gene_ID
genes_to_label <- genes_to_label %>%
  
  distinct(
    Gene_ID,
    .keep_all = TRUE
  )


# ----------------------------------------------------------
# Check that column names and gene IDs are now unique
# ----------------------------------------------------------

cat(
  "\nDuplicated EIN3 column names:",
  anyDuplicated(names(genes_to_label)),
  "\n"
)

cat(
  "Duplicated EIN3 Gene_IDs:",
  anyDuplicated(genes_to_label$Gene_ID),
  "\n"
)


stopifnot(
  anyDuplicated(names(genes_to_label)) == 0
)


stopifnot(
  anyDuplicated(genes_to_label$Gene_ID) == 0
)


cat(
  "Total unique EIN3 genes in reference list:",
  nrow(genes_to_label),
  "\n"
)


# ==========================================================
# 5. DEFINE SIGNIFICANCE THRESHOLDS
# ==========================================================

alpha <- 0.05

lfc_cut <- 1


# ==========================================================
# 6. PREPARE DESEQ2 RESULTS
# ==========================================================

df <- df %>%
  
  mutate(
    
    Log2FC =
      as.numeric(log2FoldChange),
    
    padj =
      as.numeric(padj),
    
    # Volcano y-axis
    neglog10 =
      -log10(padj),
    
    # Cap very small FDR values at 50 for plotting
    neglog10_cap =
      pmin(neglog10, 50),
    
    # Classify DEGs
    grp =
      case_when(
        
        !is.na(padj) &
          padj < alpha &
          Log2FC > lfc_cut ~ "up",
        
        !is.na(padj) &
          padj < alpha &
          Log2FC < -lfc_cut ~ "down",
        
        TRUE ~ "ns"
      )
    
  )


# ==========================================================
# 7. ADD EIN3 ANNOTATION
# ==========================================================

df <- df %>%
  
  left_join(
    
    genes_to_label,
    
    by = c(
      "gene" = "Gene_ID"
    )
    
  )


cat(
  "\nSuccessful EIN3 annotation join.\n"
)


cat(
  "Genes with EIN3 annotation:",
  sum(!is.na(df$Gene_Name)),
  "\n"
)


# ==========================================================
# 8. MAIN DIFFERENTIAL EXPRESSION RESULT
# ==========================================================

n_up <- sum(
  df$grp == "up",
  na.rm = TRUE
)

n_down <- sum(
  df$grp == "down",
  na.rm = TRUE
)

n_deg <- n_up + n_down


cat(
  "\n========================================\n"
)

cat(
  "MAIN DIFFERENTIAL EXPRESSION RESULT\n"
)

cat(
  "========================================\n"
)


cat(
  "Upregulated genes:",
  n_up,
  "\n"
)

cat(
  "Downregulated genes:",
  n_down,
  "\n"
)

cat(
  "Total significant DEGs:",
  n_deg,
  "\n"
)


# ==========================================================
# 9. SAVE ALL SIGNIFICANT DEGS
# ==========================================================

all_sig_genes <- df %>%
  
  filter(
    grp != "ns"
  ) %>%
  
  select(
    
    gene,
    
    Log2FC,
    
    padj,
    
    neglog10,
    
    grp
    
  )


write.csv(
  
  all_sig_genes,
  
  file =
    "ERF11_ethylene_vs_air_ALL_significant_DEGs.csv",
  
  row.names = FALSE
)


# ==========================================================
# CONTROL 1: KNOWN-ANSWER CONTROL
# ==========================================================
#
# Prediction:
#
# Because EIN3 is a central regulator of ethylene signaling,
# independently identified EIN3-regulated genes should be
# overrepresented among genes responding to ethylene.
#
# ==========================================================


cat(
  "\n========================================\n"
)

cat(
  "CONTROL 1: KNOWN ANSWER - EIN3\n"
)

cat(
  "========================================\n"
)


# ==========================================================
# 10. DEFINE STATISTICALLY TESTABLE UNIVERSE
# ==========================================================

analysis_universe <- df %>%
  
  filter(
    !is.na(padj)
  )


n_testable <- nrow(
  analysis_universe
)


cat(
  "Genes with non-NA adjusted P values:",
  n_testable,
  "\n"
)


# ==========================================================
# 11. EIN3 GENES PRESENT/TESTABLE
# ==========================================================

ein3_tested <- analysis_universe %>%
  
  filter(
    gene %in% genes_to_label$Gene_ID
  )


n_ein3_tested <- nrow(
  ein3_tested
)


cat(
  "Known EIN3 genes present/testable:",
  n_ein3_tested,
  "\n"
)


if (n_ein3_tested == 0) {
  
  stop(
    "No EIN3 genes matched the DESeq2 gene IDs. Check gene-ID formatting."
  )
}


# ==========================================================
# 12. EIN3 GENES RECOVERED AS DEGS
# ==========================================================

ein3_DEGs <- ein3_tested %>%
  
  filter(
    grp != "ns"
  )


n_ein3_deg <- nrow(
  ein3_DEGs
)


cat(
  "Known EIN3 genes recovered as DEGs:",
  n_ein3_deg,
  "\n"
)


# ==========================================================
# 13. RECOVERY PERCENTAGE
# ==========================================================

ein3_recovery <-
  
  100 *
  
  n_ein3_deg /
  
  n_ein3_tested


cat(
  "EIN3 recovery percentage:",
  round(
    ein3_recovery,
    1
  ),
  "%\n"
)


# ==========================================================
# 14. UP/DOWN DIRECTION OF EIN3 DEGS
# ==========================================================

cat(
  "\nDirection of recovered EIN3 genes:\n"
)


print(
  table(
    ein3_DEGs$grp
  )
)


# ==========================================================
# 15. SAVE RECOVERED EIN3 DEGS
# ==========================================================

write.csv(
  
  ein3_DEGs %>%
    
    select(
      
      gene,
      
      Gene_Name,
      
      Log2FC,
      
      padj,
      
      grp
      
    ),
  
  file =
    "Known_answer_EIN3_DEGs.csv",
  
  row.names = FALSE
)


# ==========================================================
# 16. FISHER'S EXACT ENRICHMENT TEST
# ==========================================================
#
# Question:
#
# Are known EIN3-regulated genes more common among DEGs
# than expected from all statistically testable genes?
#
# ==========================================================


analysis_universe <- analysis_universe %>%
  
  mutate(
    
    EIN3 =
      gene %in%
      genes_to_label$Gene_ID,
    
    DEG =
      grp != "ns"
    
  )


# Force FALSE/TRUE levels so table is always 2 x 2
ein3_table <- table(
  
  EIN3_regulated =
    factor(
      analysis_universe$EIN3,
      levels = c(FALSE, TRUE)
    ),
  
  DEG =
    factor(
      analysis_universe$DEG,
      levels = c(FALSE, TRUE)
    )
)


cat(
  "\nEIN3 x DEG contingency table:\n"
)


print(
  ein3_table
)


# ==========================================================
# 17. RUN FISHER TEST
# ==========================================================

ein3_fisher <- fisher.test(
  ein3_table
)


cat(
  "\nFisher's exact test:\n"
)


print(
  ein3_fisher
)


# Extract odds ratio and P value
ein3_OR <- unname(
  ein3_fisher$estimate
)


ein3_p <- ein3_fisher$p.value


cat(
  "\nEIN3 enrichment odds ratio:",
  round(
    ein3_OR,
    2
  ),
  "\n"
)


cat(
  "EIN3 enrichment P value:",
  signif(
    ein3_p,
    4
  ),
  "\n"
)


# ==========================================================
# 18. SELECT EIN3 GENES TO LABEL ON VOLCANO
# ==========================================================
#
# Only:
# 1. significant genes
# 2. EIN3-regulated
# 3. have informative gene name
# 4. Gene_Name != Gene_ID
#
# ==========================================================


label_df <- df %>%
  
  filter(
    
    !is.na(Gene_Name),
    
    Gene_Name != "",
    
    grp != "ns",
    
    Gene_Name != gene
    
  )


n_labeled <- nrow(
  label_df
)


cat(
  "\nEIN3 genes labeled on volcano:",
  n_labeled,
  "\n"
)


cat(
  "Direction of labeled EIN3 genes:\n"
)


print(
  table(
    label_df$grp
  )
)


# Validate labels
stopifnot(
  
  all(
    
    label_df$padj < alpha &
      
      abs(
        label_df$Log2FC
      ) > lfc_cut
    
  )
)


# ==========================================================
# 19. SAVE VOLCANO-LABELED EIN3 GENES
# ==========================================================

write.csv(
  
  label_df %>%
    
    select(
      
      gene,
      
      Gene_Name,
      
      Log2FC,
      
      padj,
      
      grp
      
    ),
  
  file =
    "Volcano_labeled_EIN3_genes_ERF11_ethylene_vs_air.csv",
  
  row.names = FALSE
)


# ==========================================================
# 20. SAVE ALL SIGNIFICANT EIN3 GENES
# ==========================================================
#
# Includes genes where Gene_Name == Gene_ID.
#
# ==========================================================


all_regulatory_sig_genes <- df %>%
  
  filter(
    
    !is.na(Gene_Name),
    
    grp != "ns"
    
  ) %>%
  
  select(
    
    gene,
    
    Gene_Name,
    
    Log2FC,
    
    padj,
    
    grp
    
  )


write.csv(
  
  all_regulatory_sig_genes,
  
  file =
    "ALL_significant_EIN3_genes_ERF11_ethylene_vs_air.csv",
  
  row.names = FALSE
)


# ==========================================================
# 21. CREATE KNOWN-ANSWER SUBTITLE
# ==========================================================

known_answer_subtitle <- paste0(
  
  "Known-answer control: ",
  
  n_ein3_deg,
  
  "/",
  
  n_ein3_tested,
  
  " EIN3-regulated genes recovered as DEGs (",
  
  round(
    ein3_recovery,
    1
  ),
  
  "%); OR = ",
  
  round(
    ein3_OR,
    2
  ),
  
  ", Fisher P = ",
  
  format.pval(
    
    ein3_p,
    
    digits = 2,
    
    eps = 1e-16
    
  )
)


cat(
  "\nVolcano subtitle:\n"
)

cat(
  known_answer_subtitle,
  "\n"
)


# ==========================================================
# 22. BUILD FINAL VOLCANO PLOT
# ==========================================================

p <- ggplot(
  
  df,
  
  aes(
    
    x = Log2FC,
    
    y = neglog10_cap
    
  )
  
) +
  
  
  # All genes
  geom_point(
    
    aes(
      color = grp
    ),
    
    size = 1.8,
    
    alpha = 0.85
    
  ) +
  
  
  # Colors
  scale_color_manual(
    
    values = c(
      
      down = "blue",
      
      ns = "black",
      
      up = "red"
      
    ),
    
    guide = "none"
    
  ) +
  
  
  # log2FC thresholds
  geom_vline(
    
    xintercept = c(
      -lfc_cut,
      lfc_cut
    ),
    
    linetype = "dotted"
    
  ) +
  
  
  # FDR threshold
  geom_hline(
    
    yintercept =
      -log10(alpha),
    
    linetype = "dotted"
    
  ) +
  
  
  # Number downregulated
  annotate(
    
    "text",
    
    x = -7,
    
    y = 30,
    
    label = n_down,
    
    color = "blue",
    
    fontface = "bold",
    
    size = 6
    
  ) +
  
  
  # Number upregulated
  annotate(
    
    "text",
    
    x = 7,
    
    y = 30,
    
    label = n_up,
    
    color = "red",
    
    fontface = "bold",
    
    size = 6
    
  ) +
  
  
  # EIN3-regulated gene labels
  geom_text_repel(
    
    data = label_df,
    
    aes(
      label = Gene_Name
    ),
    
    size = 3,
    
    max.overlaps = Inf,
    
    nudge_x =
      ifelse(
        
        label_df$grp == "up",
        
        2,
        
        -2
        
      ),
    
    direction = "y",
    
    segment.color =
      "grey40",
    
    segment.size =
      0.4,
    
    arrow =
      grid::arrow(
        
        length =
          grid::unit(
            
            0.015,
            
            "npc"
            
          )
        
      ),
    
    force = 2
    
  ) +
  
  
  labs(
    
    x =
      "log2(Fold Change)",
    
    y =
      "-log10(FDR)",
    
    title =
      "Ethylene ERF11 vs Air ERF11",
    
    subtitle =
      known_answer_subtitle,
    
    caption =
      "Labels represent independently identified EIN3-regulated genes."
    
  ) +
  
  
  coord_cartesian(
    
    xlim = c(
      -15,
      15
    ),
    
    ylim = c(
      0,
      50
    )
    
  ) +
  
  
  theme_classic(
    base_size = 12
  ) +
  
  
  theme(
    
    plot.title =
      element_text(
        
        size = 12,
        
        face = "bold"
        
      ),
    
    plot.subtitle =
      element_text(
        
        size = 9
        
      ),
    
    plot.caption =
      element_text(
        
        size = 8
        
      )
    
  )


# ==========================================================
# 23. DISPLAY FINAL VOLCANO
# ==========================================================

p


# ==========================================================
# 24. SAVE FINAL VOLCANO
# ==========================================================

ggsave(
  
  filename =
    "Volcano_ERF11_Ethylene_vs_Air_KnownAnswer.png",
  
  plot = p,
  
  width = 11,
  
  height = 7.5,
  
  dpi = 300
)


# ==========================================================
# 25. SAVE KNOWN-ANSWER SUMMARY
# ==========================================================

sink(
  "Known_answer_EIN3_summary.txt"
)


cat(
  "CONTROL 1: KNOWN-ANSWER CONTROL\n"
)

cat(
  "================================\n\n"
)


cat(
  "Comparison: Ethylene ERF11 vs Air ERF11\n\n"
)


cat(
  "Prediction:\n"
)

cat(
  "Independently identified EIN3-regulated genes should be\n"
)

cat(
  "overrepresented among genes responding to ethylene.\n\n"
)


cat(
  "DEG criterion:\n"
)

cat(
  "FDR <",
  alpha,
  "\n"
)

cat(
  "|log2FoldChange| >",
  lfc_cut,
  "\n\n"
)


cat(
  "Total statistically testable genes:",
  n_testable,
  "\n"
)


cat(
  "Total significant DEGs:",
  n_deg,
  "\n"
)


cat(
  "Upregulated DEGs:",
  n_up,
  "\n"
)


cat(
  "Downregulated DEGs:",
  n_down,
  "\n\n"
)


cat(
  "Unique EIN3 genes in reference list:",
  nrow(genes_to_label),
  "\n"
)


cat(
  "EIN3 genes present/testable:",
  n_ein3_tested,
  "\n"
)


cat(
  "EIN3 genes recovered as DEGs:",
  n_ein3_deg,
  "\n"
)


cat(
  "EIN3 recovery percentage:",
  round(
    ein3_recovery,
    1
  ),
  "%\n\n"
)


cat(
  "Direction of EIN3 DEGs:\n"
)


print(
  table(
    ein3_DEGs$grp
  )
)


cat(
  "\nContingency table:\n"
)


print(
  ein3_table
)


cat(
  "\nFisher's exact enrichment test:\n"
)


print(
  ein3_fisher
)


cat(
  "\nInterpretation:\n"
)


if (
  ein3_OR > 1 &&
  ein3_p < 0.05
) {
  
  cat(
    "Known EIN3-regulated genes were significantly enriched among\n"
  )
  
  cat(
    "ethylene-responsive DEGs, supporting the known-answer control.\n"
  )
  
} else if (
  ein3_OR > 1
) {
  
  cat(
    "EIN3-regulated genes showed enrichment among DEGs, but the\n"
  )
  
  cat(
    "enrichment was not statistically significant at P < 0.05.\n"
  )
  
} else {
  
  cat(
    "EIN3-regulated genes were not enriched among the identified DEGs.\n"
  )
  
}


sink()


# ==========================================================
# 26. FINAL CONSOLE SUMMARY
# ==========================================================

cat(
  "\n========================================\n"
)

cat(
  "ANALYSIS COMPLETE\n"
)

cat(
  "========================================\n"
)


cat(
  "Total DEGs:",
  n_deg,
  "\n"
)


cat(
  "Upregulated:",
  n_up,
  "\n"
)


cat(
  "Downregulated:",
  n_down,
  "\n"
)


cat(
  "EIN3 genes tested:",
  n_ein3_tested,
  "\n"
)


cat(
  "EIN3 genes recovered:",
  n_ein3_deg,
  "\n"
)


cat(
  "EIN3 recovery:",
  round(
    ein3_recovery,
    1
  ),
  "%\n"
)


cat(
  "Fisher odds ratio:",
  round(
    ein3_OR,
    2
  ),
  "\n"
)


cat(
  "Fisher P value:",
  signif(
    ein3_p,
    4
  ),
  "\n"
)


cat(
  "\nFiles created:\n"
)

cat(
  "1. ERF11_ethylene_vs_air_ALL_significant_DEGs.csv\n"
)

cat(
  "2. Known_answer_EIN3_DEGs.csv\n"
)

cat(
  "3. Volcano_labeled_EIN3_genes_ERF11_ethylene_vs_air.csv\n"
)

cat(
  "4. ALL_significant_EIN3_genes_ERF11_ethylene_vs_air.csv\n"
)

cat(
  "5. Known_answer_EIN3_summary.txt\n"
)

cat(
  "6. Volcano_ERF11_Ethylene_vs_Air_KnownAnswer.png\n"
)

