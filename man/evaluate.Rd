\name{evaluate}
\alias{evaluate}
\alias{evaluate.SCE}
\alias{evaluate.SCA}
\title{
  Evaluate SCE and SCA Model Performance
}
\description{
  Evaluate model performance for SCE or SCA models.
}
\usage{
\method{evaluate}{SCE}(object, Testing_data, Training_data, digits = 3, ...)
\method{evaluate}{SCA}(object, Testing_data, Training_data, digits = 3, ...)
}
\arguments{
  \item{object}{An SCE or SCA model object}
  \item{Testing_data}{Testing dataset}
  \item{Training_data}{Training dataset}
  \item{digits}{Number of decimal places (default: 3)}
  \item{...}{Additional arguments}
}
\value{
  Model performance metrics.
}

\seealso{
  \code{\link{SCE}}, \code{\link{SCA}}, \code{\link{predict}}
} 