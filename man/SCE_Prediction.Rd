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
  \code{\link{SCE}}, \code{\link{Model_simulation}}
}

\examples{
\dontrun{
  ## Load SCE package and the supporting packages
  library(SCE)
  library(parallel)

  ## Set number of cores to use
  options(mc.cores = 2)

  ## Load example data
  data("Streamflow_training_10var")
  data("Streamflow_testing_10var")

  ## Define predictors and predictants
  Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
  Predictants <- c("Flow")

  ## Build the SCE model
  Model <- SCE(
    Training_data = Streamflow_training_10var,
    X = Predictors,
    Y = Predictants,
    mfeature = round(0.5 * length(Predictors)),
    Nmin = 5,
    Ntree = 48,
    alpha = 0.05,
    resolution = 1000
  )

  ## Generate predictions for test data
  predictions <- SCE_Prediction(
    X_sample = Streamflow_testing_10var,
    model = Model
  )
}
} 