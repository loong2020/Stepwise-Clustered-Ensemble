\name{Model_simulation}
\alias{Model_simulation}
\alias{SCE_Prediction}
\alias{OOB_validation}

\title{
  Perform Simulations using a Trained SCE Model
}

\description{
  This function facilitates the simulation of dependent variables based on the specified independent variables utilizing a previously trained SCE model. It generates predictions for training, out-of-bag validation, and testing datasets. The function includes comprehensive input validation and handles both single and multiple predictants.
}

\usage{
Model_simulation(Testing_data, Training_data, model)

SCE_Prediction(X_sample, model)

OOB_validation(model)
}

\arguments{
  \item{Testing_data}{
    A data.frame or matrix comprising the data that will be used to test the model.
    Must contain all the predictors used in the model and must not contain missing values.
  }
  \item{Training_data}{
    A data.frame or matrix comprising the data used for training the model.
    Must contain all the predictors used in the model and must not contain missing values.
  }
  \item{model}{
    The trained SCE model object generated through the \code{SCE} function.
    This model contains the necessary information and parameters for conducting simulations.
  }
  \item{X_sample}{
    A data.frame or matrix containing the predictor variables for which predictions are to be made.
    Must contain all the predictors used in the model and must not contain missing values.
  }
}

\value{
  For \code{Model_simulation}:
  A list containing three components:
  \itemize{
    \item Training: A data.frame containing predictions for the training dataset
    \item Validation: A data.frame containing out-of-bag (OOB) predictions
    \item Testing: A data.frame containing predictions for the testing dataset
  }
  
  For \code{SCE_Prediction}:
  A matrix containing weighted ensemble predictions for each predictant.
  
  For \code{OOB_validation}:
  A data.frame containing weighted OOB predictions for each predictant.
}

\details{
  The simulation process involves the following steps:
  \enumerate{
    \item Input validation for data frames/matrices and model structure
    \item Checking for required predictors and predictants
    \item Validation of data types and missing values
    \item Generation of predictions using weighted ensemble approach
  }
  
  Input validation includes checks for:
  \enumerate{
    \item Null values in inputs
    \item Data frame/matrix structure
    \item Empty datasets
    \item Presence of required predictors and predictants
    \item Missing values in data
    \item Numeric data types
  }
  
  The prediction process:
  \enumerate{
    \item Each tree in the ensemble makes predictions
    \item Predictions are weighted by the tree's OOB performance
    \item Weighted predictions are combined to form the ensemble prediction
  }
}

\author{
  Kailong Li <lkl98509509@gmail.com>
}

\examples{
  ## Load SCE package and the supporting packages
  library(SCE)
  library(parallel)
  library(hydroGOF)

  ## Load training and testing data
  data("Training_input")
  data("Testing_input")

  ## Define independent (x) and dependent (y) variables
  Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2")
  Predictants <- c("swvl3","swvl4")

  ## Build the SCE model
  Model <- SCE(
    Training_data = Training_input,
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
    Testing_data = Testing_input,
    Training_data = Training_input,
    model = Model
  )

  ## Access the different prediction components
  training_predictions <- Results$Training
  validation_predictions <- Results$Validation
  testing_predictions <- Results$Testing
}
