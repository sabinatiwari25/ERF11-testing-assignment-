############################################################
## CONTROL 7: DETERMINISM
## ERF11 RNA-seq DESeq2 analysis
##
## Biological comparison:
## Ethylene-treated ERF11 vs Air-treated ERF11
##
## PURPOSE:
##
## Test whether the exact same input data and exact same
## DESeq2 workflow produce exactly the same result when
## independently executed twice.
##
## Expected result:
##
## Run 1 and Run 2 should have:
##
## - identical genes
## - identical baseMean
## - identical log2FoldChange
## - identical standard errors
## - identical test statistics
## - identical P values
## - identical adjusted P values
## - identical DEG classifications
##
############################################################


# ==========================================================
# 1. START FROM A CLEAN WORKSPACE
# ==========================================================
#
# WHY?
#
# We do not want this analysis to depend on objects left over
# from an earlier R session.
# ==========================================================

rm(list = ls())


# ==========================================================
# 2. SET WORKING DIRECTORY
# ==========================================================

setwd(
  "C:/Users/sabin/OneDrive - Auburn University/Computational Biology Coloquium/Assignment 2"
)


cat("\nCurrent working directory:\n")
print(getwd())


# ==========================================================
# 3. LOAD DESEQ2
# ==========================================================
#
# DESeq2 must be explicitly loaded so the analysis can run
# from a fresh R session.
# ==========================================================

library(DESeq2)


cat("\nDESeq2 version:\n")
print(packageVersion("DESeq2"))


cat("\nR version:\n")
print(R.version.string)


# ==========================================================
# 4. RECORD INPUT-FILE CHECKSUMS
# ==========================================================
#
# WHY?
#
# A checksum is a fingerprint of a file.
#
# If these checksums remain the same, we know the raw input
# files themselves have not changed between analyses.
# ==========================================================

input_files <- c(
  "air and ethylene erf11 gene count matrix.csv",
  "Metadata ERF11.txt"
)


input_md5 <- tools::md5sum(
  input_files
)


cat("\nInput-file MD5 checksums:\n")
print(input_md5)


write.table(
  input_md5,
  file = "Determinism_input_MD5.txt",
  quote = FALSE,
  col.names = FALSE
)


# ==========================================================
# 5. READ RAW COUNT DATA FROM DISK
# ==========================================================
#
# IMPORTANT:
#
# We are deliberately reading the original raw-count file.
#
# We are NOT using:
#
# - RPKM
# - previously normalized counts
# - an old countdata object
# - an old dds object
#
# This makes the test reproducible from the original input.
# ==========================================================

countdata <- as.matrix(
  read.csv(
    "air and ethylene erf11 gene count matrix.csv",
    row.names = "gene_id",
    check.names = FALSE
  )
)


cat(
  "\nGenes in raw count matrix:",
  nrow(countdata),
  "\n"
)

cat(
  "Samples in raw count matrix:",
  ncol(countdata),
  "\n"
)


# ==========================================================
# 6. CLEAN GENE IDs
# ==========================================================
#
# Use exactly the same gene-ID cleaning as the original
# DESeq2 analysis.
# ==========================================================

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
# 7. READ METADATA FROM DISK
# ==========================================================

coldata <- read.table(
  "Metadata ERF11.txt",
  header = TRUE,
  row.names = 1
)


cat("\nMetadata:\n")
print(coldata)


# ==========================================================
# 8. MATCH COUNT MATRIX TO METADATA
# ==========================================================

countdata <- countdata[
  ,
  rownames(coldata),
  drop = FALSE
]


# ==========================================================
# 9. REPEAT YOUR INPUT INVARIANTS
# ==========================================================
#
# WHY?
#
# Before testing determinism, we first make sure both runs
# are starting from valid and correctly matched inputs.
# ==========================================================

cat("\n========================================\n")
cat("INPUT CHECKS\n")
cat("========================================\n")


cat(
  "Number of samples:",
  ncol(countdata),
  "\n"
)


cat(
  "All counts non-negative:",
  all(countdata >= 0),
  "\n"
)


cat(
  "All counts integer-valued:",
  all(countdata == round(countdata)),
  "\n"
)


cat(
  "Count/metadata order identical:",
  identical(
    colnames(countdata),
    rownames(coldata)
  ),
  "\n"
)


cat(
  "Duplicated sample names:",
  anyDuplicated(
    colnames(countdata)
  ),
  "\n"
)


cat(
  "Duplicated gene IDs:",
  anyDuplicated(
    rownames(countdata)
  ),
  "\n"
)


stopifnot(
  ncol(countdata) == 6
)


stopifnot(
  all(countdata >= 0)
)


stopifnot(
  all(countdata == round(countdata))
)


stopifnot(
  identical(
    colnames(countdata),
    rownames(coldata)
  )
)


stopifnot(
  anyDuplicated(
    colnames(countdata)
  ) == 0
)


stopifnot(
  anyDuplicated(
    rownames(countdata)
  ) == 0
)


# ==========================================================
# 10. DEFINE TREATMENT ORDER
# ==========================================================
#
# Air is the reference.
#
# Therefore:
#
# positive log2FC = higher in ethylene
#
# negative log2FC = lower in ethylene
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
  sum(
    coldata$Treatment == "air_ERF11"
  ) == 3
)


stopifnot(
  sum(
    coldata$Treatment == "ethylene_ERF11"
  ) == 3
)


# ==========================================================
# 11. CREATE ONE FUNCTION CONTAINING THE COMPLETE PIPELINE
# ==========================================================
#
# WHY?
#
# If Run 1 and Run 2 were written as separate blocks of code,
# we could accidentally change a parameter between them.
#
# By putting the analysis inside one function, both runs use
# exactly the same commands.
# ==========================================================

run_ERF11_DESeq2 <- function(
    counts_input,
    metadata_input
) {
  
  
  # --------------------------------------------------------
  # Create DESeq2 dataset
  # --------------------------------------------------------
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    
    countData = round(
      counts_input
    ),
    
    colData = metadata_input,
    
    design = ~ Treatment
  )
  
  
  # --------------------------------------------------------
  # Same filtering used in original analysis
  #
  # Keep genes with:
  #
  # total counts > 20
  # --------------------------------------------------------
  
  dds <- dds[
    rowSums(
      DESeq2::counts(dds)
    ) > 20,
  ]
  
  
  # --------------------------------------------------------
  # Run DESeq2
  # --------------------------------------------------------
  
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )
  
  
  # --------------------------------------------------------
  # Explicit Ethylene vs Air contrast
  # --------------------------------------------------------
  
  res <- DESeq2::results(
    
    dds,
    
    contrast = c(
      "Treatment",
      "ethylene_ERF11",
      "air_ERF11"
    ),
    
    alpha = 0.10
  )
  
  
  # Convert to standard data frame
  result_df <- as.data.frame(
    res
  )
  
  
  # Add gene ID explicitly
  result_df$gene <- rownames(
    result_df
  )
  
  
  # Keep columns in fixed order
  result_df <- result_df[
    ,
    c(
      "gene",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    )
  ]
  
  
  # Sort genes so row order is guaranteed
  result_df <- result_df[
    order(
      result_df$gene
    ),
  ]
  
  
  rownames(result_df) <- NULL
  
  
  return(
    result_df
  )
}


# ==========================================================
# 12. RUN ANALYSIS NUMBER 1
# ==========================================================

cat("\n========================================\n")
cat("RUNNING DESEQ2 ANALYSIS 1\n")
cat("========================================\n")


run1 <- run_ERF11_DESeq2(
  countdata,
  coldata
)


cat(
  "Run 1 genes:",
  nrow(run1),
  "\n"
)


# ==========================================================
# 13. RUN ANALYSIS NUMBER 2
# ==========================================================
#
# This reconstructs the complete DESeq2 analysis again.
# ==========================================================

cat("\n========================================\n")
cat("RUNNING DESEQ2 ANALYSIS 2\n")
cat("========================================\n")


run2 <- run_ERF11_DESeq2(
  countdata,
  coldata
)


cat(
  "Run 2 genes:",
  nrow(run2),
  "\n"
)


# ==========================================================
# 14. CHECK DIMENSIONS
# ==========================================================

same_dimensions <- identical(
  dim(run1),
  dim(run2)
)


cat(
  "\nDimensions identical:",
  same_dimensions,
  "\n"
)


# ==========================================================
# 15. CHECK GENE IDs AND ORDER
# ==========================================================

same_genes <- identical(
  run1$gene,
  run2$gene
)


cat(
  "Gene IDs/order identical:",
  same_genes,
  "\n"
)


# ==========================================================
# 16. EXACT OBJECT COMPARISON
# ==========================================================
#
# identical() is very strict.
#
# TRUE means the two result data frames are exactly the same
# in their stored values and structure.
# ==========================================================

exact_identical <- identical(
  run1,
  run2
)


cat(
  "Entire result table exactly identical:",
  exact_identical,
  "\n"
)


# ==========================================================
# 17. CHECK MAXIMUM NUMERICAL DIFFERENCE
# ==========================================================
#
# Even if identical() were FALSE because of an attribute,
# we can ask whether any numerical result changed.
# ==========================================================

numeric_columns <- c(
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)


max_difference <- data.frame(
  
  Variable = numeric_columns,
  
  Maximum_absolute_difference =
    NA_real_
)


for (i in seq_along(numeric_columns)) {
  
  
  var <- numeric_columns[i]
  
  
  x <- run1[[var]]
  
  y <- run2[[var]]
  
  
  finite_values <- is.finite(x) &
    is.finite(y)
  
  
  if (any(finite_values)) {
    
    
    max_difference$Maximum_absolute_difference[i] <-
      
      max(
        abs(
          x[finite_values] -
            y[finite_values]
        )
      )
    
    
  } else {
    
    
    max_difference$Maximum_absolute_difference[i] <-
      
      NA_real_
    
  }
}


cat(
  "\nMaximum numerical differences between Run 1 and Run 2:\n"
)


print(
  max_difference
)


# ==========================================================
# 18. CHECK NA PATTERNS
# ==========================================================
#
# DESeq2 can return NA adjusted P values for some genes.
#
# Determinism requires that the same genes receive NA values
# in both runs.
# ==========================================================

same_pvalue_NA <- identical(
  
  is.na(
    run1$pvalue
  ),
  
  is.na(
    run2$pvalue
  )
)


same_padj_NA <- identical(
  
  is.na(
    run1$padj
  ),
  
  is.na(
    run2$padj
  )
)


cat(
  "\nP-value NA pattern identical:",
  same_pvalue_NA,
  "\n"
)


cat(
  "Adjusted-P NA pattern identical:",
  same_padj_NA,
  "\n"
)


# ==========================================================
# 19. DEFINE YOUR ORIGINAL DEG CRITERION
# ==========================================================
#
# Same criterion used throughout the assignment:
#
# padj < 0.05
#
# AND
#
# |log2FC| > 1
# ==========================================================

classify_DEG <- function(result) {
  
  
  ifelse(
    
    !is.na(result$padj) &
      result$padj < 0.05 &
      result$log2FoldChange > 1,
    
    "up",
    
    
    ifelse(
      
      !is.na(result$padj) &
        result$padj < 0.05 &
        result$log2FoldChange < -1,
      
      "down",
      
      "ns"
    )
  )
}


run1_class <- classify_DEG(
  run1
)


run2_class <- classify_DEG(
  run2
)


# ==========================================================
# 20. CHECK DEG CLASSIFICATIONS
# ==========================================================

same_DEG_calls <- identical(
  run1_class,
  run2_class
)


cat(
  "\nDEG classifications identical:",
  same_DEG_calls,
  "\n"
)


cat(
  "\nRun 1 DEG counts:\n"
)


print(
  table(run1_class)
)


cat(
  "\nRun 2 DEG counts:\n"
)


print(
  table(run2_class)
)


# ==========================================================
# 21. COUNT UP AND DOWN DEGs
# ==========================================================

run1_up <- sum(
  run1_class == "up"
)


run1_down <- sum(
  run1_class == "down"
)


run2_up <- sum(
  run2_class == "up"
)


run2_down <- sum(
  run2_class == "down"
)


cat(
  "\nRun 1:",
  run1_up,
  "up;",
  run1_down,
  "down;",
  run1_up + run1_down,
  "total DEGs\n"
)


cat(
  "Run 2:",
  run2_up,
  "up;",
  run2_down,
  "down;",
  run2_up + run2_down,
  "total DEGs\n"
)


# ==========================================================
# 22. SAVE BOTH RESULTS
# ==========================================================

write.csv(
  run1,
  "Determinism_DESeq2_run1.csv",
  row.names = FALSE
)


write.csv(
  run2,
  "Determinism_DESeq2_run2.csv",
  row.names = FALSE
)


# ==========================================================
# 23. COMPARE OUTPUT FILE CHECKSUMS
# ==========================================================
#
# If the CSV files have exactly the same MD5 value, they are
# byte-for-byte identical.
# ==========================================================

output_md5 <- tools::md5sum(
  c(
    "Determinism_DESeq2_run1.csv",
    "Determinism_DESeq2_run2.csv"
  )
)


cat(
  "\nOutput CSV MD5 checksums:\n"
)


print(
  output_md5
)


same_output_MD5 <-
  length(
    unique(
      unname(
        output_md5
      )
    )
  ) == 1


cat(
  "Output CSV files byte-for-byte identical:",
  same_output_MD5,
  "\n"
)


# ==========================================================
# 24. DETERMINE OVERALL PASS/FAIL
# ==========================================================

numeric_identical <- all(
  
  max_difference$Maximum_absolute_difference[
    !is.na(
      max_difference$Maximum_absolute_difference
    )
  ] == 0
)


determinism_pass <-
  
  same_dimensions &&
  
  same_genes &&
  
  same_pvalue_NA &&
  
  same_padj_NA &&
  
  same_DEG_calls &&
  
  numeric_identical


cat("\n========================================\n")


cat(
  "DETERMINISM CONTROL:",
  ifelse(
    determinism_pass,
    "PASS",
    "CHECK REQUIRED"
  ),
  "\n"
)


cat("========================================\n")


# ==========================================================
# 25. SAVE SOFTWARE ENVIRONMENT
# ==========================================================
#
# WHY?
#
# Code alone is not always enough for reproducibility.
#
# Package and R versions may affect results, so they should
# be recorded.
# ==========================================================

capture.output(
  sessionInfo(),
  file =
    "SessionInfo_DESeq2_determinism.txt"
)


# ==========================================================
# 26. SAVE PACKAGE VERSION TABLE
# ==========================================================

packages_used <- c(
  "DESeq2",
  "edgeR",
  "ggplot2",
  "dplyr",
  "ggrepel"
)


package_versions <- data.frame(
  
  Package =
    packages_used,
  
  Version =
    sapply(
      
      packages_used,
      
      function(pkg) {
        
        as.character(
          packageVersion(pkg)
        )
        
      }
    )
)


write.csv(
  package_versions,
  "Software_versions.csv",
  row.names = FALSE
)


cat(
  "\nSoftware versions:\n"
)


print(
  package_versions
)


# ==========================================================
# 27. SAVE TEXT SUMMARY
# ==========================================================

sink(
  "Determinism_summary.txt"
)


cat(
  "CONTROL 7: DETERMINISM\n"
)


cat(
  "======================\n\n"
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
  "Determine whether the same raw input data and same DESeq2\n"
)


cat(
  "pipeline produce identical differential-expression results\n"
)


cat(
  "when independently executed twice.\n\n"
)


cat(
  "Input files:\n"
)


print(
  input_md5
)


cat(
  "\nRun 1 genes:",
  nrow(run1),
  "\n"
)


cat(
  "Run 2 genes:",
  nrow(run2),
  "\n\n"
)


cat(
  "Dimensions identical:",
  same_dimensions,
  "\n"
)


cat(
  "Gene order identical:",
  same_genes,
  "\n"
)


cat(
  "Entire result table exactly identical:",
  exact_identical,
  "\n"
)


cat(
  "P-value NA pattern identical:",
  same_pvalue_NA,
  "\n"
)


cat(
  "Adjusted-P NA pattern identical:",
  same_padj_NA,
  "\n"
)


cat(
  "DEG classifications identical:",
  same_DEG_calls,
  "\n"
)


cat(
  "Output CSV files byte-for-byte identical:",
  same_output_MD5,
  "\n\n"
)


cat(
  "Maximum numerical differences:\n"
)


print(
  max_difference
)


cat(
  "\nRun 1 DEG counts:\n"
)


print(
  table(run1_class)
)


cat(
  "\nRun 2 DEG counts:\n"
)


print(
  table(run2_class)
)


cat(
  "\nOverall determinism result:",
  ifelse(
    determinism_pass,
    "PASS",
    "CHECK REQUIRED"
  ),
  "\n"
)


sink()


# ==========================================================
# 28. FINAL OUTPUT
# ==========================================================

cat("\n========================================\n")

cat("DETERMINISM ANALYSIS COMPLETE\n")

cat("========================================\n")


cat(
  "Dimensions identical:",
  same_dimensions,
  "\n"
)


cat(
  "Gene IDs/order identical:",
  same_genes,
  "\n"
)


cat(
  "Entire result table exactly identical:",
  exact_identical,
  "\n"
)


cat(
  "P-value NA pattern identical:",
  same_pvalue_NA,
  "\n"
)


cat(
  "Adjusted-P NA pattern identical:",
  same_padj_NA,
  "\n"
)


cat(
  "DEG classifications identical:",
  same_DEG_calls,
  "\n"
)


cat(
  "Numerical results identical:",
  numeric_identical,
  "\n"
)


cat(
  "Output CSV files byte-for-byte identical:",
  same_output_MD5,
  "\n"
)


cat(
  "\nRun 1 total DEGs:",
  run1_up + run1_down,
  "\n"
)


cat(
  "Run 2 total DEGs:",
  run2_up + run2_down,
  "\n"
)


cat(
  "\nDETERMINISM CONTROL:",
  ifelse(
    determinism_pass,
    "PASS",
    "CHECK REQUIRED"
  ),
  "\n"
)


cat(
  "\nFiles created:\n"
)


cat(
  "1. Determinism_DESeq2_run1.csv\n"
)


cat(
  "2. Determinism_DESeq2_run2.csv\n"
)


cat(
  "3. Determinism_input_MD5.txt\n"
)


cat(
  "4. Determinism_summary.txt\n"
)


cat(
  "5. SessionInfo_DESeq2_determinism.txt\n"
)


cat(
  "6. Software_versions.csv\n"
)