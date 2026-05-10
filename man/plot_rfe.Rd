\name{plot_rfe}
\alias{plot_rfe}
\alias{Plot_RFE}
\title{
  Plot Recursive Feature Elimination Results
}
\description{
  Plot Recursive Feature Elimination results.
}
\usage{
plot_rfe(rfe_result,
         main = NULL,
         col_validation = "blue",
         col_testing = "red",
         pch = 16,
         lwd = 2,
         cex = 1.2,
         cex_text = 1,
         legend_pos = "bottom",
         ...)
}
\arguments{
  \item{rfe_result}{Result object from \code{\link{rfe_sce}} function}
  \item{main}{Plot title. Defaults to \code{NULL} (no title) for clean publication output; pass a character string to add one.}
  \item{col_validation}{Color for validation line}
  \item{col_testing}{Color for testing line}
  \item{pch}{Point character}
  \item{lwd}{Line width}
  \item{cex}{Point size}
  \item{cex_text}{Text size used uniformly for axis tick labels, axis titles, and legend text, so the plot uses a single consistent font size.}
  \item{legend_pos}{Legend position. The default \code{"bottom"} places the legend horizontally below the plotting area so it does not overlap with the data; any other valid \code{\link[graphics]{legend}} keyword (e.g. \code{"topright"}, \code{"bottomleft"}) places it inside the panel at that position.}
  \item{...}{Additional arguments}
}
\value{
  Plot showing validation and testing \eqn{R^2} vs number of predictors.
}

\seealso{
  \code{\link{rfe_sce}}
}
