\name{SCE_Prediction}
\alias{SCE_Prediction}

\title{
  Generate Predictions using an SCE Model
}

\description{
  This function generates predictions for new data using a trained SCE model.
  It combines predictions from individual trees in the ensemble, weighted by
  their respective importance weights.
}

\usage{
SCE_Prediction(X_sample, model)
}

\arguments{
  \item{X_sample}{
    A data.frame or matrix containing the predictor variables for which
    predictions are to be made. Must contain all predictors used in model training.
  }
  \item{model}{
    A trained SCE model object generated through the \code{SCE} function.
    The model must contain the trained trees and their weights.
  }
}

\value{
  A matrix containing the ensemble predictions for each predictant.
  The number of rows equals the number of samples in X_sample, and the columns
  correspond to the predictant variables.
}

\details{
  The prediction process involves:
  \enumerate{
    \item Generating predictions from each tree in the ensemble
    \item Weighting the predictions based on tree weights
    \item Combining the weighted predictions to form the ensemble prediction
  }
  
  This function is typically called internally by \code{Model_simulation}
  and is not usually called directly by users.
}

\author{
  Kailong Li <lkl98509509@gmail.com>
}

\seealso{
  \code{\link{SCE}}
}
