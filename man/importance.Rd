\name{importance}
\alias{importance}
\alias{importance.SCA}
\alias{importance.SCE}
\alias{Wilks_importance}
\alias{SCA_importance}

\title{
  Variable Importance Analysis
}

\description{
  Functions for calculating variable importance scores using Wilks' Lambda method. The package provides both generic S3 methods and direct function calls for importance analysis.
}

\usage{
importance(object, ...)

\method{importance}{SCA}(object, ...)

\method{importance}{SCE}(object, OOB_weight = TRUE, ...)

Wilks_importance(model, OOB_weight = TRUE)

SCA_importance(model)
}

\arguments{
  \item{object}{
    An object for which importance scores should be calculated.
  }
  \item{model}{
    A trained model object:
    \itemize{
      \item For \code{Wilks_importance}: SCE model object (S3 class "SCE") containing a list of SCA trees
      \item For \code{SCA_importance}: Single SCA tree object (S3 class "SCA")
    }
  }
  \item{OOB_weight}{
    Logical indicating whether to use out-of-bag weighting for importance calculation.
    Default is TRUE. Only used for SCE objects and \code{Wilks_importance}.
  }
  \item{...}{
    Additional arguments passed to methods.
  }
}

\value{
  A data.frame containing:
  \itemize{
    \item \code{Predictor}: Names of the predictors
    \item \code{Relative_Importance}: Normalized importance scores (sum to 1)
  }
}

\details{
  \strong{Importance Calculation Method:}
  
  All functions use the Wilks' Lambda statistic to calculate variable importance:
  \enumerate{
    \item Extract Wilks' Lambda values and split information from tree(s)
    \item Replace negative Wilks' Lambda values with zero
    \item Calculate raw importance for each split:
      \itemize{
        \item Importance = (left_samples + right_samples) / total_samples * (1 - Wilks' Lambda)
      }
    \item Aggregate importance scores by predictor
    \item Normalize importance scores to sum to 1
  }
  
  \strong{Function Differences:}
  \itemize{
    \item \code{importance.SCA()}: S3 method for single SCA trees
    \item \code{importance.SCE()}: S3 method for SCE ensembles (calls \code{Wilks_importance})
    \item \code{Wilks_importance()}: Direct function for SCE ensembles with OOB weighting options
    \item \code{SCA_importance()}: Direct function for single SCA trees
  }
  
  \strong{OOB Weighting:}
  \itemize{
    \item If \code{OOB_weight = TRUE}: Importance scores are weighted by each tree's OOB performance
    \item If \code{OOB_weight = FALSE}: Importance scores are calculated using the median across trees
  }
}

\author{
  Kailong Li <lkl98509509@gmail.com>
}

\references{
  Li, Kailong, Guohe Huang, and Brian Baetz. "Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling." Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.
}

\seealso{
  \code{\link{SCA}}, \code{\link{SCE}}
}

