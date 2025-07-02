\name{Model_simulation}
\alias{Model_simulation}

\title{
  Perform Simulations using a Trained SCE Model
}

\description{
  This function facilitates the simulation of dependent variables based on the specified independent variables utilizing a previously trained SCE model. It generates predictions for training, out-of-bag validation, and testing datasets. The function includes comprehensive input validation for data types, missing values, and predictor matching.
}

\usage{
Model_simulation(model, Testing_data)
}

\arguments{
  \item{model}{
    The trained SCE model object (S3 class "SCE") generated through the \code{SCE} function.
    This model contains the necessary information and parameters for conducting simulations.
  }
  \item{Testing_data}{
    A data.frame or matrix comprising the data that will be used to test the model.
    Must contain all the predictors used in the model.
    Must not contain missing values.
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

\seealso{
  \code{\link{SCE}}
  See the generic function \code{predict} for SCE objects.
}
