# ERF11 RNA-seq Testing-Mindset Assignment

## Biological question

This analysis examines transcriptional differences between
ethylene-treated ERF11 and air-treated ERF11 samples.

## Experimental design

Six biological samples were analyzed.

### Air-treated ERF11

- air_ERF11.1
- air_ERF11.2
- air_ERF11.3

### Ethylene-treated ERF11

- ethylene_ERF11.1
- ethylene_ERF11.2
- ethylene_ERF11.3

All samples represent the ERF11 genotype.

## Input data

The primary analysis was performed using raw gene-level RNA-seq counts.

Input files:

- air and ethylene erf11 gene count matrix.csv
- Metadata ERF11.txt

The raw/unpublished RNA-seq data are not included in the public GitHub
repository.

Input-file MD5 checksums are recorded in:

Determinism_input_MD5.txt

## Primary DESeq2 analysis

Differential expression was analyzed using DESeq2.

Model:

design = ~ Treatment

Comparison:

ethylene_ERF11 vs air_ERF11

Genes with total counts <= 20 across the six samples were removed.

After filtering:

19,820 genes were analyzed.

DESeq2 results were extracted using the original independent-filtering
setting of alpha = 0.10.

Differentially expressed genes were defined using:

adjusted P value < 0.05

and

|log2FoldChange| > 1

The primary analysis identified:

- 121 upregulated genes
- 124 downregulated genes
- 245 total DEGs


# Testing Controls

## Control 1: Known-answer control

An independently defined set of EIN3-regulated genes was used as a
biological known-answer benchmark because EIN3 is a major regulator
of ethylene-responsive transcription.

Of 250 testable EIN3-regulated genes, 40 were identified as DEGs.

Fisher's exact test showed strong enrichment:

Odds ratio = 9.76

P = 2.58 × 10^-23

This supported the ability of the analysis pipeline to recover an
independently expected component of the ethylene response.


## Control 2: Invariant checks

The following invariants were tested:

- exactly six samples
- three Air and three Ethylene biological replicates
- non-negative raw counts
- integer-valued raw counts
- identical count-matrix and metadata sample order
- no duplicated sample names
- no duplicated gene IDs

All required invariant checks passed.


## Control 3: Negative control

The association between treatment labels and RNA-seq expression
profiles was deliberately disrupted while retaining the original
samples, genes, library sizes, and balanced 3-vs-3 experimental
structure.

Because six samples allow only 10 unique balanced 3-vs-3 partitions
when complementary assignments are treated as equivalent, all nine
incorrect partitions were evaluated.

Real treatment assignment:

245 DEGs

Incorrect assignments:

Median = 2 DEGs

Range = 0-83 DEGs

The real treatment assignment therefore produced substantially more
differential-expression signal than any incorrect assignment.


## Control 4: Positive-control titration

Artificial expression effects were introduced into 20 genes that
showed little evidence of differential expression in the original
analysis.

Planted effects:

- log2FC = 0.25
- log2FC = 0.5
- log2FC = 1
- log2FC = 2
- log2FC = 3

Recovery using the primary DEG criterion was:

- 0.25: 0%
- 0.5: 0%
- 1: 65%
- 2: 100%
- 3: 100%

Detection therefore increased with increasing signal strength, as
expected.


## Control 5: Redundancy

The expression response was independently estimated using edgeR.

edgeR used:

- filterByExpr
- TMM normalization
- negative-binomial dispersion estimation
- quasi-likelihood testing

15,812 genes were retained by edgeR.

Across genes tested by both methods:

Pearson correlation = 0.9849

Spearman correlation = 0.9958

Overall direction agreement = 98.4%

Among DESeq2-significant genes retained by edgeR:

Direction agreement = 100%

Pearson correlation = 0.9996

Spearman correlation = 0.9997

edgeR identified only one gene meeting the final FDR and fold-change
threshold, demonstrating that significance calls were sensitive to the
statistical framework even though effect-size estimates were highly
concordant.


## Control 6: Order-of-magnitude expectations

Expected ranges and sanity checks are documented in:

07_order_of_magnitude_expectations.md

The observed DEG counts, negative-control results, positive-control
dose response, and cross-method fold-change agreement were within
biologically and computationally plausible ranges.


## Control 7: Determinism

The complete DESeq2 workflow was executed twice independently from
the original raw count and metadata files.

Both runs produced:

- 19,820 analyzed genes
- 121 upregulated genes
- 124 downregulated genes
- 245 total DEGs

The two runs had:

- identical gene IDs and ordering
- identical baseMean values
- identical log2 fold changes
- identical standard errors
- identical test statistics
- identical P values
- identical adjusted P values
- identical DEG classifications
- identical output-file MD5 checksums

Maximum numerical difference between the two runs was zero for every
DESeq2 result variable.

Therefore, the DESeq2 workflow passed the determinism control.


# Software environment

The analysis was performed using:

- R 4.4.1
- DESeq2 1.44.0
- edgeR 4.2.2
- ggplot2 4.0.1
- dplyr 1.2.0
- ggrepel 0.9.6

Detailed software information is recorded in:

- SessionInfo_DESeq2_determinism.txt
- Software_versions.csv


# Reproducibility

Analysis scripts, testing controls, documentation, and selected outputs
are version-controlled using Git and stored in this GitHub repository.

Raw unpublished RNA-seq data are excluded from version control.