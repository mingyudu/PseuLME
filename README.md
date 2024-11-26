# Tutorial of PseuLME

## Overview

**PseuLME** is an R package designed for differential expression analysis of single-cell omics data using a pseudobulk approach combined with linear mixed-effects models. 

---

## Installation

PseuLME depends on several other R packages, including `DESeq2`, `SummarizedExperiment`, `lme4`, `multcomp`, and `stats`. After installing the dependencies, you can install the package as follows:

```r
if(!require(devtools))
  install.packages(devtools)
# Install PseuLME from GitHub
devtools::install_github("mingyudu/PseuLME")
# Load the package
library(PseuLME)
```

---

## How to use

After installing the package, you can view the documentation for the main and helper functions using:

```r
# Main function
?fit_PseuLME

# Helper functions
?create_subclusters
?filter_genes
?filter_low_celltypes
?convert_to_pseudobulk_modified
?lme_pvals
```

To conduct differential expression analysis using PseuLME, prepare a `DESeqDataSet` object (`dds`) containing your raw count matrix at single cell level and corresponding metadata. Then, run the algorithm. Here is an example code:

```r
data(dds)
res = fit_PseuLME(dds, 
                  cell_type_accessor = 'cluster_coarse', 
                  condition_accessor = 'condition', 
                  sample_accessor = 'sample_id', 
                  batch_accessor = 'assay_id', 
                  gene_pct = 0.05, 
                  sample_thres = 50, 
                  contrast1 = 'Ir-Veh', 
                  contrast2 = 'Ir-Rz')
```

---

## Parameter Explanation

- **`dds`**: A `DESeqDataSet` object containing raw count data for all cells and associated metadata. The count matrix should be provided as rows (genes) and columns (cells), with metadata columns specifying attributes like cell type, condition, sample, and batch, etc.

- **`cell_type_accessor`**: The name of the column in the metadata (`colData(dds)`) that identifies cell types. This column is used to create subclusters for differential expression analysis (e.g., `"cluster_coarse"`).

- **`condition_accessor`**: The name of the column in the metadata (`colData(dds)`) that identifies experimental conditions, such as control and treatment groups (e.g., `"condition"`).

- **`sample_accessor`**: The name of the column in the metadata (`colData(dds)`) that specifies sample identifiers. This is used for pseudobulk aggregation (e.g., `"sample_id"`).

- **`batch_accessor`**: The name of the column in the metadata (`colData(dds)`) that specifies batch identifiers. This accounts for batch effects in the linear mixed-effects model (e.g., `"assay_id"`).

- **`gene_pct`**: A numeric threshold (between 0 and 1) representing the minimum proportion of cells in which a gene must be expressed to be retained. For example, `gene_pct = 0.05` keeps genes expressed in at least 5% of cells.

- **`sample_thres`**: An integer threshold specifying the minimum number of cells required in a sample to include it in the analysis. For example, `sample_thres = 50` retains samples with at least 50 cells.

- **`contrast1`**: The label of the control group in the condition column (e.g., `"Ir-Veh"`).

- **`contrast2`**: The label of the case group in the condition column (e.g., `"Ir-Rz"`).

---

## Output

The result, `res`, is a named list where each element corresponds to a cell type subcluster. Each element is a data frame with the following columns for each gene:

- **`base`**: Baseline expression level.
- **`log2fc`**: Log2 fold change between `contrast2` and `contrast1`.
- **`pvals`**: P-values from the differential expression test.
- **`padjs`**: FDR-adjusted p-values using the Benjamini-Hochberg method.
- **`gene`**: Name of genes.
