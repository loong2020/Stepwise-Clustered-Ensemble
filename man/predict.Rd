\name{predict}
\alias{predict}
\alias{predict.sce}
\alias{predict.sca}
\title{
  Predict Using SCE and SCA Models
}
\description{
  Make predictions on new data using SCE or SCA models.
}
\usage{
\method{predict}{sce}(object, newdata, ...)
\method{predict}{sca}(object, newdata, ...)
}
\arguments{
  \item{object}{An SCE or SCA model object}
  \item{newdata}{New data for prediction}
  \item{...}{Additional arguments}
}
\value{
  Predictions for the new data.
}

\seealso{
  \code{\link{sce}}, \code{\link{sca}}, \code{\link{evaluate}}
} 