\name{print.SCA}
\alias{print.SCA}
\alias{summary.SCA}
\alias{print.SCE}
\alias{summary.SCE}

\title{
  Print and Summary Methods for SCA and SCE Objects
}

\description{
  Methods for printing and summarizing SCA (Stepwise Cluster Analysis) and SCE (Stepwise Clustered Ensemble) objects.
}

\usage{
\method{print}{SCA}(x, ...)

\method{summary}{SCA}(object, ...)

\method{print}{SCE}(x, ...)

\method{summary}{SCE}(object, ...)
}

\arguments{
  \item{x, object}{
    An SCA or SCE object returned by the \code{SCA()} or \code{SCE()} function.
  }
  \item{...}{
    Additional arguments passed to methods.
  }
}

\value{
  \code{print.SCA} and \code{print.SCE} return the object invisibly.
  \code{summary.SCA} returns a summary of the SCA model including tree statistics.
  \code{summary.SCE} returns a summary of the SCE model including ensemble statistics.
}

\seealso{
  \code{\link{SCA}}, \code{\link{SCE}}
}
 