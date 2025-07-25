\name{predict}
\alias{predict}
\alias{predict.SCE}
\alias{predict.SCA}
\title{
  Predict Using SCE and SCA Models
}
\description{
  Make predictions on new data using SCE or SCA models.
}
\usage{
\method{predict}{SCE}(object, newdata, ...)
\method{predict}{SCA}(object, newdata, ...)
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
  \code{\link{SCE}}, \code{\link{SCA}}, \code{\link{evaluate}}
} 