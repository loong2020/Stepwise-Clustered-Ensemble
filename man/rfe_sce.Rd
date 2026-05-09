\name{rfe_sce}
\alias{rfe_sce}
\alias{RFE_SCE}
\title{
  Recursive Feature Elimination for SCE Models
}
\description{
  Recursive Feature Elimination for SCE models to identify the most important predictors.
}
\usage{
rfe_sce(training_data, testing_data, predictors, predictant, nmin, ntree,
        alpha = 0.05, resolution = 1000, step = 1, verbose = TRUE,
        parallel = TRUE)
}
\arguments{
  \item{training_data}{Training dataset}
  \item{testing_data}{Testing dataset}
  \item{predictors}{Character vector of predictor names}
  \item{predictant}{Character vector of predictant names}
  \item{nmin}{Minimum samples per node}
  \item{ntree}{Number of trees}
  \item{alpha}{Significance level (default: 0.05)}
  \item{resolution}{Resolution for splitting (default: 1000)}
  \item{step}{Number of predictors to remove per iteration (default: 1)}
  \item{verbose}{Print progress (default: TRUE)}
  \item{parallel}{Use parallel processing (default: TRUE)}
}
\value{
  RFE results with performance metrics and importance scores.
}

\seealso{
  \code{\link{plot_rfe}}, \code{\link{sce}}, \code{\link{importance}}
}
