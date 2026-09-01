############################################################
## Volcano plot of DESeq2 results
## Comparison: air wrky estradiol vs air wrky
## Purpose:
##  - Plot all genes from DESeq2
##  - Highlight significant up/down regulated genes
##  - Cap -log10(FDR) at 50 (flat top line)
##  - Label ONLY significant regulatory genes
##  - Do NOT label genes where Gene_Name == Gene_ID
############################################################

# ----------------------------------------------------------
# 1. Set working directory (where all files are stored)
# ----------------------------------------------------------
setwd("C:/Users/sabin/OneDrive - Auburn University/Computational Biology Coloquium/Assignment 2")

# ----------------------------------------------------------
# 2. Load required libraries
# ----------------------------------------------------------
library(ggplot2)   # plotting
library(dplyr)     # data wrangling
library(ggrepel)   # non-overlapping text labels


# ----------------------------------------------------------
# 3. Read input files
# ----------------------------------------------------------

# Main DESeq2 result file
df <- read.csv("1SEPDGESeq_results_ethylene_ERF11_vs_air_ERF11.csv")

# Regulatory gene list (must contain Gene_ID and Gene_Name columns)
genes_to_label <- read.csv("EIN3_ethylene.csv")


# ----------------------------------------------------------
# 4. Define significance thresholds
# ----------------------------------------------------------
alpha   <- 0.05   # FDR cutoff
lfc_cut <- 1      # log2 fold change cutoff


# ----------------------------------------------------------
# 5. Prepare DESeq2 data for plotting
# ----------------------------------------------------------
df <- df %>%
  mutate(
    # Gene ID column (original DESeq2 IDs)
    gene = X,
    
    # Ensure numeric values
    Log2FC = as.numeric(log2FoldChange),
    padj   = as.numeric(padj),
    
    # Convert FDR to -log10 scale for volcano plot y-axis
    neglog10 = -log10(padj),
    
    # Cap values above 50 so very significant genes
    # align on a straight horizontal line
    neglog10_cap = pmin(neglog10, 50),
    
    # Classify genes into up / down / non-significant
    grp = case_when(
      padj < alpha & Log2FC >  lfc_cut ~ "up",
      padj < alpha & Log2FC < -lfc_cut ~ "down",
      TRUE                             ~ "ns"
    )
  ) %>%
  # Join regulatory gene annotation
  left_join(genes_to_label, by = c("gene" = "Gene_ID"))


# ----------------------------------------------------------
# 6. Count number of significant genes
# ----------------------------------------------------------
n_up   <- sum(df$grp == "up",   na.rm = TRUE)
n_down <- sum(df$grp == "down", na.rm = TRUE)

# Create table of ALL significant DEGs (up + down)
# ----------------------------------------------------------

library(dplyr)     # data wrangling
all_sig_genes <- df %>%
  dplyr::filter(grp != "ns") %>%
  dplyr::select(
    gene,
    Log2FC,
    padj,
    neglog10,
    grp
  )

# Save to CSV
write.csv(
  all_sig_genes,
  file = "ethyl ERF11 vs air ERF11_ALL_sig_genes.csv",
  row.names = FALSE
)

# ----------------------------------------------------------
# 7. Select genes to label on the plot
# ----------------------------------------------------------
label_df <- df %>%
  filter(
    !is.na(Gene_Name),   # gene is in regulatory list
    grp != "ns",         # gene is significant (up or down)
    Gene_Name != gene    # exclude cases where name == ID
  )

label_df %>%
  select(
    gene,
    Gene_Name,
    Log2FC,
    padj,
    grp
  ) %>%
  write.csv(
    file = "Volcano_labeled_genes_ethyl ERF11 vs air ERF11.csv",
    row.names = FALSE
  )
#label_df <- df %>%
#filter(
# !is.na(Gene_Name),
# grp != "ns",
#Gene_Name != gene
#) %>%
#arrange(padj) %>%     # most significant first
#slice(1:20)           # label top 20 only

# ----------------------------------------------------------
# 8. Build the volcano plot
# ----------------------------------------------------------
p <- ggplot(df, aes(x = Log2FC, y = neglog10_cap)) +
  
  # Plot all genes as points
  geom_point(
    aes(color = grp),
    size  = 1.8,
    alpha = 0.85
  ) +
  
  # Color scheme for groups
  scale_color_manual(
    values = c(
      down = "blue",
      ns   = "black",
      up   = "red"
    ),
    guide = "none"
  ) +
  
  # Vertical lines for fold-change cutoffs
  geom_vline(
    xintercept = c(-lfc_cut, lfc_cut),
    linetype = "dotted"
  ) +
  
  # Horizontal line for FDR cutoff
  geom_hline(
    yintercept = -log10(alpha),
    linetype = "dotted"
  ) +
  
  # Annotate number of down-regulated genes
  annotate(
    "text",
    x = -7, y = 30,
    label = n_down,
    color = "blue",
    fontface = "bold",
    size = 6
  ) +
  
  # Annotate number of up-regulated genes
  annotate(
    "text",
    x = 7, y = 30,
    label = n_up,
    color = "red",
    fontface = "bold",
    size = 6
  ) +
  
  # Label significant regulatory genes only
  geom_label_repel(
    data = label_df,
    aes(label = Gene_Name),
    size          = 3,
    max.overlaps  = Inf,
    box.padding   = 0.4,                      # space around label box
    point.padding = 0.2,                      # space around point
    label.size    = 0.2,                      # border thickness of box
    label.r       = unit(0.15, "lines"),      # roundness of box corners
    segment.color = "grey40",
    segment.size  = 0.4,
    arrow         = arrow(length = unit(0.015, "npc"))  # <- arrow to point
  ) +
  
  labs(
    x = "log2(Fold Change)",
    y = "-log10(FDR)",
    title = "ethylene ERF11 vs air ERF11"
  ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 50)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 10, face = "bold"))



# ----------------------------------------------------------
# 9. Display the plot
# ----------------------------------------------------------
p

# Number of genes labeled on the volcano plot
n_labeled <- nrow(label_df)
n_labeled

# Check that all labeled genes are significant (FDR < 0.05 and |log2FC| > 1)
all(label_df$padj < 0.05 & abs(label_df$Log2FC) > 1)

# Count labeled genes that are up- vs down-regulated
table(label_df$grp)

p <- ggplot(df, aes(x = Log2FC, y = neglog10_cap)) +
  geom_point(aes(color = grp), size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(down = "blue", ns = "black", up = "red"), guide = "none") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
  geom_hline(yintercept = -log10(alpha), linetype = "dotted") +
  annotate("text", x = -7, y = 30, label = n_down, color = "blue", fontface = "bold", size = 6) +
  annotate("text", x = 7, y = 30, label = n_up, color = "red", fontface = "bold", size = 6) +
  geom_text_repel(
    data = label_df,
    aes(label = Gene_Name),
    size = 3,
    max.overlaps = Inf,
    nudge_x = ifelse(label_df$grp == "up", 2.0, -2.0),
    direction = "y",
    segment.color = "grey40",
    segment.size  = 0.4,
    arrow = arrow(length = unit(0.015, "npc")),
    force = 2
  ) +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(FDR)",
    title = "ethylene ERF11 vs air ERF11"
  ) +
  coord_cartesian(xlim = c(-15, 15), ylim = c(0, 50)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 10, face = "bold"))

p

# Save ALL significant regulatory genes
# Includes genes where Gene_Name == Gene_ID
all_regulatory_sig_genes <- df %>% filter( !is.na(Gene_Name), # gene is in regulatory list 
                                           grp != "ns" # significant up- or downregulated gene
) %>% select( gene, Gene_Name, Log2FC, padj, grp ) 
write.csv( all_regulatory_sig_genes,
           file = "ALL_significant_EIN3_ethylene_genes_ethylene ERF11 vs air ERF11.csv", row.names = FALSE )
