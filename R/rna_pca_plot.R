#' @title PCA plot of the RNA samples
#'
#'@description rna_pca_plot function makes a pca plot using the PC1 and PC2 to see how the samples are distributed.
#'
#' @param pca_dat data frame that is the result after performing pca_prep function
#' @param x PC that should be displayed on the x axis. Default is 1, which is PC1.
#' @param y PC that should be displayed on the y axis. Default is 2, which is PC2.
#'
#' @return
#' @export
#'
#' @examples
#' pca_data_nb <- pca_prep(norm_dat = sample_vst_dd, "Group")
#'
#' # To plot PC1 and PC2
#' ana <- "Sample Title PC1 vs PC2"
#' caption <- "Data Source: \n Sample Data"
#'
#' rna_pca_plot(pca_dat = pca_data_nb, x = 1, y = 2)
#'
#' # To plot PC2 and PC3
#' rna_pca_plot(pca_dat = pca_data_nb, x = 2, y = 3)
#'
rna_pca_plot <- function(pca_dat, x = 1, y = 2) {

  if (x == 1 & y == 2) {
    xlab <- paste0("PC1: ",pca_dat$percentVar[1],"% variance")
    ylab <- paste0("PC2: ",pca_dat$percentVar[2],"% variance")

    p <- ggplot(pca_dat, aes(x=PC1, y=PC2, color = condition))

  }

  if (x == 2 & y == 3) {
    xlab <- paste0("PC2: ",pca_dat$percentVar[2],"% variance")
    ylab <- paste0("PC3: ",pca_dat$percentVar[3],"% variance")

    p <- ggplot(pca_dat, aes(x=PC2, y=PC3, color = condition))

  }

  if (x == 3 & y == 4) {
    xlab <- paste0("PC3: ",pca_dat$percentVar[3],"% variance")
    ylab <- paste0("PC4: ",pca_dat$percentVar[4],"% variance")

    p <- ggplot(pca_dat, aes(x=PC3, y=PC4, color = condition))

  }

  if (x == 4 & y == 5) {
    xlab <- paste0("PC4: ",pca_dat$percentVar[4],"% variance")
    ylab <- paste0("PC5: ",pca_dat$percentVar[5],"% variance")

    p <- ggplot(pca_dat, aes(x=PC4, y=PC5, color = condition))

  }

  if (x == 5 & y == 6) {
    xlab <- paste0("PC5: ",pca_dat$percentVar[5],"% variance")
    ylab <- paste0("PC6: ",pca_dat$percentVar[6],"% variance")

    p <- ggplot(pca_dat, aes(x=PC5, y=PC6, color = condition))

  }

  p + geom_point(size = 2.2, alpha = 0.7) +
    ggrepel::geom_text_repel(data=pca_dat, aes(label=sampleNO), size = 3.5)+
    xlab(xlab) + ylab(ylab) +
    theme_bw() +
    scale_fill_brewer(type="qual") +
    scale_colour_brewer(type="qual", palette=7) +
    scale_x_continuous(breaks=seq(-100,100,5)) +
    scale_y_continuous(breaks=seq(-100,100,5)) +
    labs(title = ana, caption = caption,
         subtitle = "PCs calculated from top 1000 variable genes.")

}
