############################################################
## CONTROL 3: NEGATIVE CONTROL
## ERF11 Air vs Ethylene RNA-seq
##
## Destroy:
##   True relationship between sample and treatment
##
## Preserve:
##   Raw counts
##   Gene identities
##   Library sizes
##   Six samples
##   3-vs-3 group sizes
##   DESeq2 pipeline
##   DEG thresholds
############################################################

setwd(
  "C:/Users/sabin/OneDrive - Auburn University/Computational Biology Coloquium/Assignment 2"
)

library(DESeq2)
library(dplyr)
library(ggplot2)
packageVersion("DESeq2")
exists("DESeqDataSetFromMatrix")

countdata <- as.matrix(
  read.csv(
    "air and ethylene erf11 gene count matrix.csv",
    row.names = "gene_id",
    check.names = FALSE
  )
)

# Clean gene IDs exactly as in original analysis
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

coldata <- read.table(
  "Metadata ERF11.txt",
  header = TRUE,
  row.names = 1
)

countdata <- countdata[, rownames(coldata)]
identical(
  colnames(countdata),
  rownames(coldata)
)
real_treatment <- as.character(
  coldata$Treatment
)

real_treatment

#identify the sample names:
samples <- colnames(countdata)

samples

#Generate all possible 3-sample groupings
all_groups <- combn(
  samples,
  3,
  simplify = FALSE
)

length(all_groups)

#keep only one member of each complementary pair
unique_groups <- all_groups[
  sapply(
    all_groups,
    function(x) samples[1] %in% x
  )
]

length(unique_groups)

#Identify the real biological partition
real_air <- rownames(coldata)[
  coldata$Treatment == "air_ERF11"
]

real_air

#determine which of the 10 partitions corresponds to the actual Air/Ethylene experiment:
is_real <- sapply(
  unique_groups,
  function(x) {
    setequal(x, real_air) ||
      setequal(setdiff(samples, x), real_air)
  }
)

which(is_real)

#create only the intentionally incorrect assignments:
null_groups <- unique_groups[!is_real]

length(null_groups)

#same filtering and thresholds as your real analysis:
alpha <- 0.05
lfc_cut <- 1
negative_results <- data.frame(
  permutation = integer(),
  n_up = integer(),
  n_down = integer(),
  n_DEG = integer()
)

#run all nine negative controls:
for (i in seq_along(null_groups)) {
  
  group1 <- null_groups[[i]]
  
  perm_meta <- data.frame(
    Treatment = factor(
      ifelse(
        samples %in% group1,
        "air_ERF11",
        "ethylene_ERF11"
      ),
      levels = c(
        "air_ERF11",
        "ethylene_ERF11"
      )
    ),
    row.names = samples
  )
  
  # Create DESeq2 dataset
  dds_perm <- DESeqDataSetFromMatrix(
    countData = round(countdata),
    colData = perm_meta,
    design = ~ Treatment
  )
  
  # SAME filtering rule as real analysis
  dds_perm <- dds_perm[
    rowSums(counts(dds_perm)) > 20,
  ]
  
  # Run DESeq2 quietly
  dds_perm <- DESeq(
    dds_perm,
    quiet = TRUE
  )
  
  res_perm <- results(
    dds_perm,
    contrast = c(
      "Treatment",
      "ethylene_ERF11",
      "air_ERF11"
    ),
    alpha = alpha
  )
  
  res_perm <- as.data.frame(
    res_perm
  )
  
  # Same DEG definition
  up <- sum(
    !is.na(res_perm$padj) &
      res_perm$padj < alpha &
      res_perm$log2FoldChange > lfc_cut
  )
  
  down <- sum(
    !is.na(res_perm$padj) &
      res_perm$padj < alpha &
      res_perm$log2FoldChange < -lfc_cut
  )
  
  negative_results <- rbind(
    negative_results,
    data.frame(
      permutation = i,
      n_up = up,
      n_down = down,
      n_DEG = up + down
    )
  )
  
  cat(
    "Negative control",
    i,
    ":",
    up + down,
    "DEGs\n"
  )
}

#Look at your null distribution
negative_results
summary(
  negative_results$n_DEG
)
real_DEGs <- 245
real_DEGs
max(negative_results$n_DEG)
median(negative_results$n_DEG)

#Save the negative-control results
write.csv(
  negative_results,
  "Negative_control_permutations.csv",
  row.names = FALSE
)

#Make the negative-control figure
p_neg <- ggplot(
  negative_results,
  aes(x = factor(permutation), y = n_DEG)
) +
  
  geom_col() +
  
  geom_hline(
    yintercept = real_DEGs,
    linetype = "dashed"
  ) +
  
  annotate(
    "text",
    x = 5,
    y = real_DEGs,
    label = paste0(
      "Real treatment = ",
      real_DEGs,
      " DEGs"
    ),
    vjust = -0.2
  ) +
  
  labs(
    x = "Incorrect treatment assignment",
    y = "Number of significant DEGs",
    title = "Negative-control treatment-label permutations",
    subtitle =
      "True Air/Ethylene relationship intentionally destroyed"
  ) +
  
  theme_classic(base_size = 13)

p_neg

ggsave(
  "Negative_control_permutations.png",
  p_neg,
  width = 8,
  height = 6,
  dpi = 300
)

negative_results

summary(negative_results$n_DEG)

median(negative_results$n_DEG)

max(negative_results$n_DEG)
