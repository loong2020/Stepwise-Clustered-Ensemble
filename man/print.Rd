\name{print}
\alias{print}
\alias{print.sce}
\alias{print.sca}
\title{
  Print SCE and SCA Model Objects
}
\description{
  Print information about SCE or SCA model objects.
}
\usage{
\method{print}{sce}(x, ...)
\method{print}{sca}(x, ...)
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
  \code{\link{sce}}, \code{\link{sca}}, \code{\link{summary}}
} 