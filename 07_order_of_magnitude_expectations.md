# ERF11 RNA-seq Testing Assignment

## Biological comparison

Differential gene expression between ethylene-treated ERF11
and air-treated ERF11 samples.

## Samples

Six biological samples were analyzed:

Ethylene:
- ethylene_ERF11.1
- ethylene_ERF11.2
- ethylene_ERF11.3

Air:
- air_ERF11.1
- air_ERF11.2
- air_ERF11.3

All samples represent the ERF11 genotype.

## Primary input

Raw gene-level RNA-seq count matrix:

air and ethylene erf11 gene count matrix.csv

Metadata:

Metadata ERF11.txt

Raw/unpublished data are excluded from the public Git repository.

## Primary differential-expression analysis

DESeq2 was used with:

design = ~ Treatment

Comparison:

ethylene_ERF11 versus air_ERF11

Genes with total counts <= 20 across the six samples were removed.

DEG criterion:

FDR < 0.05
and
|log2FoldChange| > 1

## Testing controls

1. Known-answer control using independently identified EIN3-regulated genes.
2. Input/data invariants.
3. Exhaustive treatment-label negative control.
4. Artificial positive-control signal titration.
5. Redundancy analysis using edgeR.
6. Order-of-magnitude sanity checks.
7. Deterministic rerunning of the DESeq2 workflow.

## Reproducibility

Software versions are recorded in:

SessionInfo_DESeq2_determinism.txt

and

Software_versions.csv


