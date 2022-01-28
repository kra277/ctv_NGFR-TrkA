#' Calculate inflation factor
#'
#' @description The function calculates the inflation factor.
#' This function is part of the Q-Q plot generation and is needed by the \code{\link{rna_qq_plot}}.
#'
#' @param pvals Unadjusted P values
#'
#' @return
#' @export
#'
#' @examples
#' inflation(deg_df$pvalue)
#'
inflation <- function(pvals) {
  chisq <- qchisq(1 - pvals, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
