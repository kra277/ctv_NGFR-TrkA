#' Quantile - Quantile plots for RNA seq
#'
#' @description Function for creating the Quantile - Quantile plots of P values from DEG data frames
#'
#' @param pvals Un adjusted p values to make the QQ plot
#' @param ci Confidence interval. Default is 0.95.
#'
#' @return
#' @export
#'
#' @examples
#' caption <- "Data Source: \n Sample Data"
#' rna_qq_plot(sample_clean_deg$pvalue)
#'
#'
rna_qq_plot <- function(pvals, ci = 0.95) {
  n <- length(pvals)
  df <- data.frame(
    observed = -log10(sort(pvals)),
    expected = -log10(ppoints(n)),
    clower = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw(base_size = 15) +
    labs(caption = caption,
         subtitle = "Observed unadjusted P values are plotted against the expected.") +
    annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.15,
      vjust = 1 + 0.15 * 3,
      label = sprintf("Inflation Factor = %.2f", inflation(pvals)),
      size = 6.2) +
    theme(
      axis.ticks = element_line(size = 0.5),
      panel.grid = element_blank()
    )
}
