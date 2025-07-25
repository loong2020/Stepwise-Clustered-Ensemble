\name{print}
\alias{print}
\alias{print.SCE}
\alias{print.SCA}
\title{
  Print SCE and SCA Model Objects
}
\description{
  Print information about SCE or SCA model objects.
}
\usage{
\method{print}{SCE}(x, ...)
\method{print}{SCA}(x, ...)
}
\arguments{
  \item{x}{An SCE or SCA model object}
  \item{...}{Additional arguments (not used)}
}
\value{
  Prints model information and returns the object invisibly.
}
\details{
  For SCE objects, prints ensemble information including number of trees, parameters, predictors, predictants, and OOB performance metrics.

  For SCA objects, prints tree structure information including total nodes, leaf nodes, cutting/merging actions, and variable names.
}
\seealso{
  \code{\link{SCE}}, \code{\link{SCA}}, \code{\link{summary}}
} 