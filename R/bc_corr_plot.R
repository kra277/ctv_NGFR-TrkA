
#' Correlation plots of the surrogate variables
#'
#'
#'@description Function for generating correlation plots for the surrogate variables and the variable of interest.
#'
#' @param bc_method Which batch correction method was performed. Two option, ruv and sva.
#' @param var_int is the variable of interest, this is assuming you only have one variable of interest.
#' @param font_size correlation font size inside each box.
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' # After running SVA
#'
#'
#' bc_corr_plot(bc_method = "sva", var_int = "Group", font_size = 1)
#'
#'
bc_corr_plot <-
  function(dat, bc_method, var_int = var_int, font_size) {

    # Get all the Phenotypes of interest
    pheno <-
      colData %>%
      dplyr::rename("variable"= {{var_int}}) %>%
      dplyr::select(variable) %>%
      dplyr::mutate(variable = as.numeric(variable)) %>%
      drop_na()


    if (bc_method == "ruv") {

      ### Prep the data

      # Add the surrogate variables as a dataframe
      ruv_df <- as.data.frame(set@phenoData@data)

      # Rename all the columns
      names(ruv_df) <- gsub("W_", "RUV_SV", names(ruv_df))

      # Merge the sv_df with he phenotype
      pheno_ruv <- cbind(pheno, ruv_df)

      # Generating the correlation matrix
      corr_ruv <<- cor(pheno_ruv)


      col1 <- colorRampPalette(c("#2B226D", "#372C8C", "#4235A9",
                                 "#6259AF", "#B5B1D3",
                                 "#D7B49E", "#DC602E",
                                 "#D4592E", "#CC512D", "#BC412B"))

      # Generate correlation plot
      corrplot(corr_ruv, method = "number",
               type = "full",
               title = "Correlation between the RUV surrogate variables and phenotype",
               mar = c(0, 0, 1, 1),
               col = col1(10),
               number.cex = font_size,
               number.digits = 2,
               tl.col = "black")

    }

    if (bc_method == "sva") {

      # Add the surrogate variables as a dataframe
      sv_df <- as.data.frame(svseq$sv)

      # Rename all the columns
      names(sv_df) <- gsub("V", "SV", names(sv_df))

      sv_df <<- sv_df

      # Merge the sv_df with he phenotype
      pheno_sv <- cbind(pheno, sv_df)

      # Generating the correlation matrix
      corr_sva <<- cor(pheno_sv)

      col1 <- colorRampPalette(c("#2B226D", "#372C8C", "#4235A9",
                                 "#6259AF", "#B5B1D3",
                                 "#D7B49E", "#DC602E",
                                 "#D4592E", "#CC512D", "#BC412B"))

      # Generate correlation plot
      corrplot(corr_sva, method = "number",
               type = "full",
               title = "Correlation between the SVA Surrogate variables and Phenotypes",
               mar = c(0, 0, 1, 1),
               col = col1(10),
               number.cex = 0.7,
               number.digits = 2,
               tl.col = "black")

    }

}

