\name{SCA}
\alias{SCA}
\alias{SCA_tree_predict}

\title{
  Stepwise Cluster Analysis (SCA) Model
}

\description{
  This function implements a Stepwise Cluster Analysis (SCA) model for multivariate data analysis. The SCA model recursively partitions the data space based on Wilks' Lambda statistic, creating a tree structure that can be used for prediction. The function includes comprehensive input validation for data types, missing values, and sample size requirements, and supports both single and multiple predictants.
}

\usage{
SCA(Training_data, X, Y, Nmin, alpha = 0.05, resolution = 1000)

SCA_tree_predict(Test_data, X, model)
}

\arguments{
  \item{Training_data}{
    A data.frame or matrix containing the training data.
    Must include all specified predictors and predictants.
    Must not contain missing values.
  }
  \item{X}{
    A character vector specifying the names of independent variables
    (e.g., c("Prcp","SRad","Tmax")). Must be present in Training_data.
    All variables must be numeric.
  }
  \item{Y}{
    A character vector specifying the name(s) of dependent variable(s)
    (e.g., c("swvl3","swvl4")). Must be present in Training_data.
    All variables must be numeric.
  }
  \item{Nmin}{
    Integer specifying the minimal number of samples in a leaf node for cutting.
    Must be a positive number and less than the sample size.
  }
  \item{alpha}{
    Numeric significance level for clustering, between 0 and 1.
    Default value is 0.05.
  }
  \item{resolution}{
    Numeric value specifying the resolution for splitting.
    Controls the granularity of the search for optimal split points.
    Default value is 1000.
  }
  \item{Test_data}{
    A data.frame or matrix containing the test data.
    Must include all predictors used in model training.
    Must not contain missing values.
  }
  \item{model}{
    A trained SCA model object returned by the SCA function.
  }
}

\value{
  For \code{SCA}:
  A list containing:
  \itemize{
    \item Tree: The SCA tree structure
    \item Map: Mapping information for predictions
    \item XName: Names of predictors used
    \item YName: Names of predictants
    \item type: Mapping type (currently "mean")
    \item totalNodes: Total number of nodes in the tree
    \item leafNodes: Number of leaf nodes
    \item cuttingActions: Number of cutting actions performed
    \item mergingActions: Number of merging actions performed
  }
  
  For \code{SCA_tree_predict}:
  A list containing:
  \itemize{
    \item All components from the input model
    \item Testing_sim: A data.frame containing predictions for the test data
  }
}

\details{
  The SCA model building process involves:
  \enumerate{
    \item Input validation:
      \itemize{
        \item Data type and structure checks (data.frame or matrix)
        \item Missing value checks in both predictors and predictants
        \item Numeric data validation for all variables
        \item Sample size requirements (must be greater than Nmin)
        \item Parameter validation (alpha, Nmin, resolution)
      }
    \item Data preparation:
      \itemize{
        \item Conversion of input data to matrix format
        \item Dimension checks and storage
        \item Parameter initialization
      }
    \item Tree construction:
      \itemize{
        \item Recursive partitioning based on Wilks' Lambda
        \item Node splitting and merging
        \item Leaf node creation
      }
  }
  
  The prediction process:
  \enumerate{
    \item Input validation for test data:
      \itemize{
        \item Data type and structure checks
        \item Missing value checks
        \item Predictor matching with training data
        \item Numeric data validation
      }
    \item Data preparation:
      \itemize{
        \item Conversion to matrix format
        \item Initialization of prediction matrix
      }
    \item Tree traversal and prediction:
      \itemize{
        \item Processing each test sample through the tree
        \item Generating predictions using leaf node mappings
      }
  }
}

\author{
  Xiuquan Wang <xxwang@upei.ca> (original SCA)
  Kailong Li <lkl98509509@gmail.com> (Resolution-search-based SCA)
}

\examples{
  ## Load SCE package
  library(SCE)

  ## Load training and testing data
  data("Training_input")
  data("Testing_input")

  ## Define independent (x) and dependent (y) variables
  Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2")
  Predictants <- c("swvl3","swvl4")

  ## Build the SCA model
  Model <- SCA(
    Training_data = Training_input,
    X = Predictors,
    Y = Predictants,
    Nmin = 5,
    alpha = 0.05,
    resolution = 1000
  )

  ## Make predictions
  Predictions <- SCA_tree_predict(
    Test_data = Testing_input,
    X = Predictors,
    model = Model
  )

  ## Access predictions
  testing_predictions <- Predictions$Testing_sim
  print("Predictions:")
  print(testing_predictions)
} 