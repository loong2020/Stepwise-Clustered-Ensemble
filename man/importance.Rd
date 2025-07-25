\name{importance}
\alias{importance}
\alias{importance.SCE}
\alias{importance.SCA}
\title{
  Variable Importance for SCE and SCA Models
}
\description{
  Calculate variable importance for SCE or SCA models.
}
\usage{
\method{importance}{SCE}(object, OOB_weight = TRUE, ...)
\method{importance}{SCA}(object, ...)
}
\arguments{
  \item{object}{An SCE or SCA model object}
  \item{OOB_weight}{Use out-of-bag weights for importance calculation (SCE only, default: TRUE)}
  \item{...}{Additional arguments}
}
\value{
  Variable importance rankings.
}

\seealso{
  \code{\link{SCE}}, \code{\link{SCA}}, \code{\link{RFE_SCE}}
} 