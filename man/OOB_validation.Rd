\name{OOB_validation}
\alias{OOB_validation}

\title{
  Calculate Out-of-Bag (OOB) Validation Predictions
}

\description{
  This function calculates the out-of-bag (OOB) validation predictions for an SCE model.
  OOB predictions are made using only the trees that did not use a particular observation
  during training, providing an unbiased estimate of the model's performance.
}

\usage{
OOB_validation(model)
}

\arguments{
  \item{model}{
    A trained SCE model object generated through the \code{SCE} function.
    The model must contain OOB information in its Tree_Info component.
  }
}

\value{
  A data.frame containing the OOB predictions for each predictant.
  The number of rows equals the number of training samples, and the columns
  correspond to the predictant variables.
}

\details{
  The OOB validation process involves:
  \enumerate{
    \item Extracting OOB indices and predictions from each tree
    \item Weighting the predictions based on tree weights
    \item Calculating weighted means for each sample
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
}
} 