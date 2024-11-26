# Tutorial of PseuLME

## Installation
PseuLME depends on a few other R packages: DESeq2, lme4, multcomp, stats, and SummarizedExperiment.
```r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("mingyudu/PseuLME")
library(PseuLME)
```

## How to use

After installing the package, check the document of the functions in this R package by running the following code:
```r
# main function
?fit_PseuLME

# other functions
?create_subclusters
?filter_genes
?filter_low_celltypes
?convert_to_pseudobulk_modified
?lme_pvals
```

After put the raw count matrix and its corresponding metadata into a DESeqDataSet `dds`, conduct the PseuLME by running the following code:
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
