#' Annotate arrange Deseq2 results
#'
#' @description Function for adding q-values (FDR), annotate, arrange, and clean the deseq2 results.
#'
#' Steps involved:
#' \enumerate{
#'
#'   \item Deseq2 results were converted into a data frame.
#'   \item Row names were converted to the column (ENTREZID).
#'   \item Q values (FDR) were calculated using the P-values and 'qvalue' package.
#'   \item ENTREZID were Annotated to get the SYMBOL and GENENAME.
#'   \item Annotations were merged with results.
#'   \item The columns were rearranged and the results were arranged based on P values.
#'
#'}
#'
#' @param deseq_res a Large Deseq Results
#'
#' @return clean_res_df: Deseq2 results which is DEG arranged ascendingly based on their P values
#'
#' @export
#'
#' @examples
#' clean_res <- clean_deseq_res(deseq_res)
#'
clean_deseq_res <- function(deseq_res) {
  ## To add q-values, annotate, arrange, and clean the deseq2 results

  deg_res_dat <-
    as.data.frame(deseq_res) %>% ## convert into a dataframe
    tibble::rownames_to_column("ENTREZID") %>% ## convert rownames to column
    dplyr::mutate(qvalue = qvalue::qvalue(p=deseq_res$pvalue)$qvalues) ## Get qvalues

  library(Homo.sapiens)

  ## Get annotations
  anno <- AnnotationDbi::select(Homo.sapiens, keytype='ENTREZID',
                                keys=deg_res_dat$ENTREZID,
                                columns=c('SYMBOL', 'GENENAME'),
                                multiVals="first")
  ## Clean and arrange
  clean_res_df <-
    deg_res_dat %>%
    merge(anno) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(ENTREZID, SYMBOL, log2FoldChange,
                  pvalue, padj, qvalue, everything())

  ## Return the clean DEG dataframe
  return(clean_res_df)
}
