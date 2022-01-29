#' @title Prep for PCA plot
#'
#' @description pca_prep function uses the vst/log transformed data to make a data frame that could be used to make a pca plot.
#'
#' @param norm_dat data frame of normalized gene counts. Generally a made using the \code{\link{DESeq2}}.
#' @param ntop number of top genes that you need to calculate the PCs. Default is 1000 genes.
#'
#' @return
#' @export
#'
#' @examples
#' pca_data_nb <- pca_prep(norm_dat = sample_vst_dd, "Group")
#'
pca_prep <- function(norm_dat, var_int, ntop = 1000) {
  # no.of genes showing the highest variance are used for PCA
  ntop = ntop

  # Rename variable of interest into var_int
  df_req <- SummarizedExperiment::colData(norm_dat) %>%
    as.data.frame() %>%
    dplyr::rename("var_int"= all_of(var_int))

  Pvars <- matrixStats::rowVars(SummarizedExperiment::assay(norm_dat))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                        length(Pvars)))]

  PCA <- prcomp(t(SummarizedExperiment::assay(norm_dat)[select, ]), scale = F)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

  # Copy the PCA Data into a Dataframe
  dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                        PC3 = PCA$x[,3], PC4 = PCA$x[,4],
                        PC5 = PCA$x[,5], PC6 = PCA$x[,6],
                        percentVar = percentVar,
                        sampleNO = colnames(norm_dat),
                        condition = df_req$var_int)

  return(dataGG)
}
