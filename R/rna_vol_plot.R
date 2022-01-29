#' @title Volcano plot for RNA DEG
#'
#' @description Function for creating the volcano plots from DEG RNA seq data frames.
#' Note: Default parameters would highlight the lo2foldchanges > |0.5|.
#'
#' @param df data frame with differential gene expression with P values and Log fold changes. Deseq2 DEG results.
#' @param p pvalue (un adjusted) or padj (adjusted p value). Default is adjusted p value (padj).
#' @param val threshold to be colored. Default is 0.1.
#'
#' @return
#' @export
#'
#' @examples
#' # for using adjusted p values and highlight values less than 0.05
#' caption <- "Data Source: \n Sample Data"
#' rna_vol_plot(clean_deg, p = padj, val = 0.05)
#'
#' # for using un - adjusted p values and highlight values less than 0.05
#' rna_vol_plot(sample_clean_deg, p = pvalue, val = 0.05)
#'
rna_vol_plot <- function(df, p = padj, val = 0.1) {

  #Colour Palette
  pal <- c(
    "Down-regulated" = "red",
    "Not significant" = "gray",
    "Up-regulated" = "blue",
    "~ Down-regulated" = "#E45093",
    "~ Up-regulated" = "#01A5E1"
  )

  # Creating thresholds
  df %>% mutate(key = factor(case_when({{p}} < val & log2FoldChange > 1 ~ "Up-regulated",
                                       {{p}} < val & log2FoldChange > 0.5 ~ "~ Up-regulated",
                                       {{p}} < val & log2FoldChange < -1 ~ "Down-regulated",
                                       {{p}} < val & log2FoldChange < -0.5 ~ "~ Down-regulated",
                                       TRUE ~ "Not significant"))) %>%
    ggplot2::ggplot(., aes(x = log2FoldChange, y = -log10({{p}}),
                  name = SYMBOL)) +
    geom_point(aes(color = key), size =1.5, alpha = 0.8)+
    scale_color_manual(values = pal)+
    geom_vline(xintercept = 0.5, color="#b8b8b8", linetype=2)+
    geom_vline(xintercept = -0.5, color="#b8b8b8", linetype=2)+
    geom_hline(yintercept = -log10(val), color="#b8b8b8", linetype=2)+
    scale_x_continuous(breaks=seq(-10,10,0.2)) +
    scale_y_continuous(breaks=seq(0,20,1)) +
    theme_bw() +
    labs(subtitle = "-log10 Un-adj P values are represented on the Y axis, and Log 2 Fold change on the X axis.",
         caption = caption)

}
