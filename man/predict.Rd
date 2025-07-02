\name{predict.SCA}
\alias{predict.SCA}
\alias{predict.SCE}
\alias{Model_simulation}
\alias{SCA_tree_predict}
\alias{SCE_Prediction}
\alias{OOB_validation}

\title{
  Model Prediction and Simulation
}

\description{
  Functions for making predictions and performing simulations using trained SCA and SCE models. The package provides both S3 methods and direct function calls for various prediction scenarios.
}

\usage{
\method{predict}{SCA}(object, newdata, ...)

\method{predict}{SCE}(object, newdata, ...)

Model_simulation(model, Testing_data)

SCA_tree_predict(model, Testing_data)

SCE_Prediction(X_sample, model)

OOB_validation(model)
}

\arguments{
  \item{object}{
    An object for which predictions should be made.
  }
  \item{newdata}{
    A data.frame or matrix containing new data for prediction.
    Must contain the same predictor variables as used in training.
  }
  \item{model}{
    A trained model object:
    \itemize{
      \item For \code{Model_simulation}: SCE model object (S3 class "SCE")
      \item For \code{SCA_tree_predict}: SCA model object (S3 class "SCA")
      \item For \code{SCE_Prediction}: SCE model object (S3 class "SCE")
      \item For \code{OOB_validation}: SCE model object (S3 class "SCE")
    }
  }
  \item{Testing_data}{
    A data.frame or matrix comprising the data that will be used to test the model.
    Must contain all the predictors used in the model.
    Must not contain missing values.
  }
  \item{X_sample}{
    A data.frame or matrix containing the predictor variables for which
    predictions are to be made. Must contain all predictors used in model training.
  }
  \item{...}{
    Additional arguments passed to methods.
  }
}

\value{
  \strong{For S3 methods:}
  \itemize{
    \item \code{predict.SCA()}: A matrix of predicted values for the predictant variables
    \item \code{predict.SCE()}: A list containing Training, Validation, and Testing predictions
  }
  
  \strong{For direct functions:}
  \itemize{
    \item \code{Model_simulation()}: A list containing three components:
      \itemize{
        \item Training: Predictions for the training dataset
        \item Validation: Out-of-bag (OOB) predictions
        \item Testing: Predictions for the testing dataset
      }
    \item \code{SCA_tree_predict()}: A list containing predictions for the test data
    \item \code{SCE_Prediction()}: A matrix containing ensemble predictions for each predictant
    \item \code{OOB_validation()}: A data.frame containing OOB predictions for each predictant
  }
}

\details{
  \strong{Prediction Methods:}
  
  \itemize{
    \item \code{predict.SCA()}: S3 method for single SCA trees (calls \code{SCA_tree_predict})
    \item \code{predict.SCE()}: S3 method for SCE ensembles (calls \code{Model_simulation})
    \item \code{Model_simulation()}: Comprehensive simulation for SCE models with training, validation, and testing predictions
    \item \code{SCA_tree_predict()}: Direct function for single SCA tree predictions
    \item \code{SCE_Prediction()}: Direct function for SCE ensemble predictions
    \item \code{OOB_validation()}: Internal function for calculating out-of-bag predictions
  }
  
  \strong{Prediction Process:}
  
  \emph{For SCA models:}
  \enumerate{
    \item Input validation (data types, missing values, predictor matching)
    \item Data preparation (conversion to matrix format)
    \item Tree traversal and prediction using leaf node mappings
  }
  
  \emph{For SCE models:}
  \enumerate{
    \item Input validation (data types, missing values, predictor matching)
    \item Data preparation (conversion to matrix format)
    \item Training predictions using all trees
    \item Out-of-bag predictions using trees not trained on each sample
    \item Testing predictions using all trees
    \item Weighting predictions based on tree weights
  }
  
  \strong{Out-of-Bag (OOB) Validation:}
  \itemize{
    \item OOB predictions are made using only trees that did not use a particular observation during training
    \item Provides unbiased estimate of model performance
    \item Used internally by \code{Model_simulation} for validation predictions
  }
  
  \strong{Input Validation:}
  All functions perform comprehensive validation:
  \enumerate{
    \item Data type and structure checks (data.frame or matrix)
    \item Missing value checks
    \item Predictor matching with training data
    \item Numeric data validation
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

