#' Create Subclusters by Cell Type
#'
#' This function splits a DESeqDataSet object into subclusters based on a
#' specified cell type accessor. Each subcluster contains data for a single cell
#' type, facilitating downstream analysis specific to each cell type.
#'
#' @param dds A `DESeqDataSet` object of data for all cell types.
#' @param accessor A character string specifying the column name in `colData(dds)`
#' that identifies cell types.
#'
#' @return A named list of `DESeqDataSet` objects, where each element corresponds
#' to a single cell type and contains only the cells of that type.
#' @export
#' @importFrom SummarizedExperiment colData
#'
create_subclusters <- function(dds,accessor) {
  subtypes=unique(colData(dds)[accessor][,1])
  subclusters=vector('list',length(subtypes))
  for (k in 1:length(subtypes)){
    # Subset the cells for the current cell type
    ind=which(colData(dds)[accessor][,1]==subtypes[k])

    subclusters[[k]]<-dds[,ind]
  }
  names(subclusters)<-subtypes
  return(subclusters)
}

#' Filter Genes by Expression Threshold for an Individual Cell Type
#'
#' This function filters genes in a DESeqDataSet based on the percentage of cells
#' in which they are expressed. It retains only those genes whose expression exceeds
#' a specified threshold, facilitating downstream analyses with more relevant genes.
#'
#' @param dds A `DESeqDataSet` object representing cells from an individual cell type.
#' @param pct A numeric value (between 0 and 1) specifying the minimum percentage of
#'    cells in which a gene must be expressed to be retained. Default is 15%.
#'
#' @return A filtered `DESeqDataSet` object containing only genes that meet the
#'    specified expression threshold.
#' @export
#' @importFrom SummarizedExperiment assays
filter_genes<-function(dds,pct=0.15){
  # Calculate the proportion of cells in which each gene is expressed
  ind=which((rowSums(assays(dds)$counts>0)/dim(dds)[2])>pct)

  dds<-dds[ind,]
  return(dds)
}

#' Filter Samples Based on Cell Counts for an Individual Cell Type
#'
#' This function filters samples in a DESeqDataSet based on the number of cells
#' for a specified cell type. It retains only samples with cell counts exceeding
#' a given threshold, facilitating analyses focused on well-represented samples.
#'
#' @param dds A `DESeqDataSet` object representing cells from an individual cell type.
#' @param sample A character string specifying the column name in `colData(dds)`
#'    that contains sample labels.
#' @param threshold An integer specifying the minimum number of cells required
#'  for a sample to be retained. Default is 50.
#'
#' @return A filtered `DESeqDataSet` object containing only samples with
#'    sufficient cell counts.
#' @export
#' @importFrom SummarizedExperiment colData
filter_low_celltypes<-function(dds,sample,threshold=50){
  good_ind=c()
  samples=unique(colData(dds)[sample][,1])
  for (k in 1:length(samples)){
    ind=which(colData(dds)[sample][,1]==samples[k])

    # Check if the cell number exceeds the threshold
    if (length(ind)>=threshold){
      good_ind=c(good_ind,ind)
    }
  }
  dds<-dds[,good_ind]
}

#' Convert Single-cell Count Matrix to Pseudobulk Matrix for a Individual Cell Type
#'
#' This function aggregates single-cell count data into pseudobulk count data
#' at the sample level for an individual cell type. It generates a `DESeqDataSet`
#' object at the pseudobulk level, which can be used for downstream differential
#' expression analyses.
#'
#' @param dds A `DESeqDataSet` object containing single-cell count data for an individual cell type.
#' @param sample_accessor A character string specifying the column name in `colData(dds)`
#'    that contains sample labels.
#' @param condition_accessor A character string specifying the column name in `colData(dds)` that identifies condition labels (e.g., treatment vs. control).
#' @param design A formula or matrix specifying the experimental design.
#'    The formula defines how the counts for each gene depend on the variables
#'    in `colData(dds)`. Default is `~1` (intercept only).
#'
#' @return A `DESeqDataSet` object with pseudobulk-level counts aggregated by sample.
#' @export
#' @import DESeq2 SummarizedExperiment
convert_to_pseudobulk_modified<-function(dds,sample_accessor,condition_accessor, design=NULL){
  dds$combined=colData(dds)[sample_accessor][,1]
  colnam=colnames(colData(dds))
  metadata=data.frame(c())

  samples=names(table(dds$combined))

  ind=which(dds$combined==samples[1])
  pseudobulk=data.frame(rowSums(assays(dds)$counts[,ind]))

  # Construct metadata
  batch_data=c()
  for (j in 1:length(colnam)){
    if (is.factor(colData(dds)[colnam[j]][,1])){
      batch_data=c(batch_data,levels(colData(dds)[colnam[j]][,1])[colData(dds)[colnam[j]][ind[1],1]])
    } else{
      batch_data=c(batch_data,colData(dds)[colnam[j]][ind[1],1])
    }
  }
  metadata=rbind(metadata,t(data.frame(batch_data)))

  # Construct pseudobulk counts
  for (j in 2:length(samples)){
    ind=which(dds$combined==samples[j])
    pseudobulk=cbind(pseudobulk,data.frame(rowSums(assays(dds)$counts[,ind])))
    batch_data=c()
    for (k in 1:length(colnam)){
      if (is.factor(colData(dds)[colnam[k]][,1])){
        batch_data=c(batch_data,levels(colData(dds)[colnam[k]][,1])[colData(dds)[colnam[k]][ind[1],1]])
      } else{
        batch_data=c(batch_data,colData(dds)[colnam[k]][ind[1],1])
      }
    }
    metadata=rbind(metadata,t(data.frame(batch_data)))
  }
  rownames(metadata)<-samples
  colnames(pseudobulk)<-samples
  colnames(metadata)<-colnam

  if(!condition_accessor %in% colnames(colData(dds))){
    stop('The condition_accessor is not in the metadata.')
  }
  if (is.null(design)){
    design=as.formula(paste0('~',condition_accessor))
  }
  dds_new<-DESeqDataSetFromMatrix(countData = pseudobulk,
                                  colData = metadata,
                                  design= design)
  colnames(dds_new)<-samples
  return(dds_new)
}

#' Fit Linear Mixed-Effect Model on Normalized Gene Expression Data
#'
#' @param dds A `DESeqDataSet` object containing pseudobulk count data for an individual cell type.
#' @param condition_accessor A character string specifying the column name in
#'    `colData(dds)` that contains condition labels (e.g., treatment vs. control).
#' @param batch_accessor A character string specifying the column name in
#'    `colData(dds)` that contains batch labels (e.g., assay 1,2,3,...).
#' @param contrast1 A character string representing the label of the control group.
#' @param contrast2 A character string representing the label of the case group.
#'
#' @return A data frame containing the following columns:
#'   - `base`: The baseline expression for each gene.
#'   - `log2fc`: The log2 fold change between `contrast2` and `contrast1`.
#'   - `pvals`: The p-value for the differential expression analysis.
#'   - `padjs`: The FDR-adjusted p-value (using the Benjamini-Hochberg method).
#'   The rows of the data frame correspond to genes in the `DESeqDataSet`.
#' @export
#' @import DESeq2 lme4 multcomp stats SummarizedExperiment
lme_pvals<-function(dds,condition_accessor,batch_accessor,
                    contrast1,contrast2){
  dds<-estimateSizeFactors(dds)
  base = rep(NA,dim(dds)[1])
  diff = rep(NA,dim(dds)[1])
  pvals=rep(NA,dim(dds)[1])
  for (k in 1:dim(dds)[1]){
    # Normalize counts using DESeq2 size factor
    data=data.frame(assays(dds)$counts[rownames(dds)[k],])/sizeFactors(dds)
    colnames(data)<-'gene'
    data$gene = log2(data$gene+1)

    # Add condition and batch information
    data$condition=factor(colData(dds)[condition_accessor][,1])
    data$batch=factor(colData(dds)[batch_accessor][,1])

    # Filter data to include only 2 contrasts
    data<-data[data$condition %in% c(contrast1,contrast2),]
    data$condition = factor(data$condition, levels = c(contrast1, contrast2))

    # model formula
    design = 'gene~1+condition+(1|batch)'

    # LME model
    temp<- tryCatch({
      lmer(formula=design,data=data)
    },error = function(e) {
      return(0)
    })
    try(base[[k]]<-summary(temp)$coefficients[1,1])
    try(post.hoc<-glht(temp,linfct=mcp(condition='Tukey')))
    try(test<-summary(post.hoc,test = adjusted("none")))

    if (!is.numeric(temp)){
      diff[[k]]<-test$test$coefficients[paste0(contrast1,' - ',contrast2)]
      pvals[[k]]<-test$test$pvalues[paste0(contrast1,' - ',contrast2)]
      if (is.na(pvals[[k]])){
        diff[[k]]<-test$test$coefficients[paste0(contrast2,' - ',contrast1)]
        pvals[[k]]<-test$test$pvalues[paste0(contrast2,' - ',contrast1)]
      }

    }
  }
  pvals[is.na(pvals)]=1
  # BH multiple testing correction
  padjs = p.adjust(pvals, method = 'BH')
  df = data.frame(base = base, log2fc = diff, pvals = pvals, padjs = padjs)
  rownames(df) = rownames(dds)
  df$gene = rownames(df)
  df = df[order(df$padjs),]
  return(df)
}

#' Cell Type Specific Pseudobulk Differential Expression Analysis via Linear Mixed-Effects Model
#'
#' This function performs differential expression analysis for subclusters of single-cell data by:
#' 1. Creating subclusters based on a specified cell type.
#' 2. Filtering genes and samples based on user-defined thresholds.
#' 3. Aggregating single-cell data into pseudobulk counts for each sample.
#' 4. Fitting linear mixed-effects models to estimate differential expression between two conditions.
#'
#' @param dds A `DESeqDataSet` object containing single-cell or pseudobulk count data.
#' @param cell_type_accessor A character string specifying the column name in `colData(dds)` that identifies cell types.
#' @param condition_accessor A character string specifying the column name in `colData(dds)` that identifies condition labels (e.g., treatment vs. control).
#' @param sample_accessor A character string specifying the column name in `colData(dds)` that identifies sample labels.
#' @param batch_accessor A character string specifying the column name in `colData(dds)` that identifies batch labels.
#' @param gene_pct A numeric value (between 0 and 1) indicating the minimum percentage of cells in which a gene must be expressed to be retained.
#' @param sample_thres An integer specifying the minimum number of cells required for a sample to be retained in the analysis.
#' @param contrast1 A character string representing the label of the control group in the condition column.
#' @param contrast2 A character string representing the label of the case group in the condition column.
#'
#' @return A named list of data frames, where each element corresponds to a cell type subcluster.
#'    Each data frame contains the following columns for each gene:
#'   - `base`: The baseline expression.
#'   - `log2fc`: The log2 fold change between `contrast2` and `contrast1`.
#'   - `pvals`: The p-value for the differential expression test.
#'   - `padjs`: The FDR-adjusted p-value (using the Benjamini-Hochberg method).
#' @export
#' @import DESeq2 lme4 multcomp stats
fit_PseuLME <- function(dds, cell_type_accessor, condition_accessor,
                   sample_accessor, batch_accessor,
                   gene_pct, sample_thres,
                   contrast1, contrast2){
  subclusters<-create_subclusters(dds, cell_type_accessor)

  # Pre-processing
  for (k in 1:length(subclusters)){
    subclusters[[k]]<-filter_genes(subclusters[[k]], gene_pct)
    subclusters[[k]]<-filter_low_celltypes(subclusters[[k]], sample_accessor, sample_thres)
  }

  # Convert to pseudobulk counts
  subclusters_pseudo<-subclusters
  rm_lst = NULL
  for (k in 1:length(subclusters)){
    try(subclusters_pseudo[[k]]<-convert_to_pseudobulk_modified(subclusters[[k]], sample_accessor, condition_accessor))
    print(names(subclusters_pseudo)[k])
    print(dim(subclusters_pseudo[[k]])[2])
    if(dim(subclusters_pseudo[[k]])[2]<6){
      rm_lst = c(rm_lst, names(subclusters_pseudo)[k])
      print(paste0(names(subclusters_pseudo)[k], ' cluster does not have enough cells, so it should be remove!'))
    }
  }

  subclusters_pseudo <- subclusters_pseudo[sapply(subclusters_pseudo, function(dds) dim(dds)[2] >= 6)]

  # DE
  res <- lapply(subclusters_pseudo, function(subcluster) {
    tryCatch(
      lme_pvals(subcluster, condition_accessor, batch_accessor, contrast1, contrast2),
      error = function(e) NULL
    )
  })
  names(res) <- names(subclusters_pseudo)
  return(res)
}


