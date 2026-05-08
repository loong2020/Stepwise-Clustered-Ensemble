\name{evaluate}
\alias{evaluate}
\alias{evaluate.sce}
\alias{evaluate.sca}
\title{
  Evaluate SCE and SCA Model Performance
}
\description{
  Evaluate model performance for SCE or SCA models.
}
\usage{
\method{evaluate}{sce}(object, Testing_data, Training_data, digits = 3, ...)
\method{evaluate}{sca}(object, Testing_data, Training_data, digits = 3, ...)
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
  \code{\link{sce}}, \code{\link{sca}}, \code{\link{predict}}
} 