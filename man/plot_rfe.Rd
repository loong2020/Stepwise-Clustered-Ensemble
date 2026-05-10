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
         legend_pos = "bottomleft",
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
  \item{legend_pos}{Legend position}
  \item{...}{Additional arguments}
}
\value{
  Plot showing validation and testing \eqn{R^2} vs number of predictors.
}

\seealso{
  \code{\link{rfe_sce}}
}
