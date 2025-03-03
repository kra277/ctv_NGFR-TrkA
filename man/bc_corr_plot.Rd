% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bc_corr_plot.R
\name{bc_corr_plot}
\alias{bc_corr_plot}
\title{Correlation plots of the surrogate variables}
\usage{
bc_corr_plot(dat, bc_method, var_int = var_int, font_size)
}
\arguments{
\item{bc_method}{Which batch correction method was performed. Two option, ruv and sva.}

\item{var_int}{is the variable of interest, this is assuming you only have one variable of interest.}

\item{font_size}{correlation font size inside each box.}
}
\value{
A correlation plot in the plottin window along with the matrix used to create the plot.
}
\description{
Function for generating correlation plots for the surrogate variables and the variable of interest.
}
\examples{

# Only execute after running RUV batch correction

\dontrun{
bc_corr_plot(bc_method = "ruv", var_int = "Group", font_size = 1)
}


# Only execute after running SVA batch correction

\dontrun{
bc_corr_plot(bc_method = "sva", var_int = "Group", font_size = 1)
}


}
