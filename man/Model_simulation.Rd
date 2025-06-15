\name{Model_simulation}
\alias{Model_simulation}

\title{
  Perform Simulations using a Trained SCE Model
}

\description{
  This function facilitates the simulation of dependent variables based on the specified independent variables utilizing a previously trained SCE model. It generates predictions for training, out-of-bag validation, and testing datasets. The function includes comprehensive input validation for data types, missing values, and predictor matching.
}

\usage{
Model_simulation(Testing_data, model)
}

\arguments{
  \item{Testing_data}{
    A data.frame or matrix comprising the data that will be used to test the model.
    Must contain all the predictors used in the model.
    Must not contain missing values.
  }
  \item{model}{
    The trained SCE model object generated through the \code{SCE} function.
    This model contains the necessary information and parameters for conducting simulations.
  }
}

\value{
  A list containing three components:
  \itemize{
    \item Training: A data.frame containing predictions for the training dataset
    \item Validation: A data.frame containing out-of-bag (OOB) predictions
    \item Testing: A data.frame containing predictions for the testing dataset
  }
}

\details{
  The simulation process involves the following steps:
  \enumerate{
    \item Input validation:
      \itemize{
        \item Data type and structure checks (data.frame or matrix)
        \item Missing value checks
        \item Predictor matching with training data
        \item Numeric data validation
      }
    \item Data preparation:
      \itemize{
        \item Conversion to matrix format
        \item Initialization of prediction matrices
      }
    \item Prediction generation:
      \itemize{
        \item Training predictions using all trees
        \item Out-of-bag predictions using trees not trained on each sample
        \item Testing predictions using all trees
      }
  }
}

\author{
  Kailong Li <lkl98509509@gmail.com>
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

  ## Define independent (x) and dependent (y) variables
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

  ## Conduct simulations using the trained SCE model
  Results <- Model_simulation(
    Testing_data = Streamflow_testing_10var,
    model = Model
  )

  ## Access the different prediction components
  training_predictions <- Results$Training
  validation_predictions <- Results$Validation
  testing_predictions <- Results$Testing
}
}
