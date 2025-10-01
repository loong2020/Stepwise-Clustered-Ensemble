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
\method{importance}{SCE}(object, OOB_weight = TRUE, digits = 2, ...)
\method{importance}{SCA}(object, digits = 2, ...)
}
\arguments{
  \item{object}{An SCE or SCA model object}
  \item{OOB_weight}{Use out-of-bag weights for importance calculation (SCE only, default: TRUE)}
  \item{digits}{Number of decimal places to round the returned relative importance values (default: 2)}
  \item{...}{Additional arguments}
}
\value{
  Variable importance rankings. For convenience, relative importance values are rounded to \code{digits} decimal places.
}

\seealso{
  \code{\link{SCE}}, \code{\link{SCA}}, \code{\link{RFE_SCE}}
} 