\name{summary}
\alias{summary}
\alias{summary.SCE}
\alias{summary.SCA}
\title{
  Summary methods for SCE and SCA models
}
\description{
  Provide concise summaries of model structure and performance for SCE and SCA objects.
}
\usage{
\method{summary}{SCE}(object, ...)
\method{summary}{SCA}(object, ...)
}
\arguments{
  \item{object}{An SCE or SCA model object}
  \item{...}{Additional arguments passed to or ignored by methods}
}
\details{
For \code{summary.SCE}, the method prints ensemble configuration, out-of-bag (OOB) performance statistics, tree structure information, and tree weight distribution.\\
For \code{summary.SCA}, the method prints tree structure information and variable summaries for the single SCA tree.
}
\value{
  Invisibly returns the input \code{object} after printing the summary.
}
\seealso{
  \code{\link{SCE}}, \code{\link{SCA}}, \code{\link{print}}, \code{\link{importance}}, \code{\link{evaluate}}
}
