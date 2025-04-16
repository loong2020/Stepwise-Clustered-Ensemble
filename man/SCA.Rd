\name{SCA}
\alias{SCA}
\alias{SCA_tree_predict}

\title{
  Stepwise Cluster Analysis (SCA) Model
}

\description{
  This function implements a Stepwise Cluster Analysis (SCA) model for multivariate data analysis. The SCA model recursively partitions the data space based on Wilks' Lambda statistic, creating a tree structure that can be used for prediction. The function includes comprehensive input validation and supports both single and multiple predictants.
}

\usage{
SCA(alpha, Training_data, X, Y, Nmin, resolution)

SCA_tree_predict(Test_data, X, model)
}

\arguments{
  \item{alpha}{
    Numeric significance level for clustering, usually between 0.001 and 0.1.
    This parameter controls the stopping criterion for tree growth.
  }
  \item{Training_data}{
    A data.frame or matrix containing the training data.
    Must include all specified predictors and predictants.
  }
  \item{X}{
    A character vector specifying the names of independent variables
    (e.g., c("Prcp","SRad","Tmax")). Must be present in Training_data.
  }
  \item{Y}{
    A character vector specifying the name(s) of dependent variable(s)
    (e.g., c("swvl3","swvl4")). Must be present in Training_data.
  }
  \item{Nmin}{
    Integer specifying the minimal number of samples in a leaf node for cutting.
    Must be greater than the number of predictants.
  }
  \item{resolution}{
    Numeric value specifying the resolution for splitting (instead of searching every possible split point, only search every resolution-th point, and iterate until the potential splits is less than resolution).
    Controls the granularity of the search for optimal split points.
  }
  \item{Test_data}{
    A data.frame or matrix containing the test data.
    Must include all predictors used in model training.
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
        \item Data type and structure checks
        \item Missing value checks
        \item Numeric data validation
        \item Sample size requirements
      }
    \item Data preparation:
      \itemize{
        \item Conversion to matrix format
        \item Dimension checks
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
    \item Input validation for test data
    \item Data preparation and initialization
    \item Tree traversal for each test sample
    \item Prediction using leaf node mappings
  }
  
  The function handles:
  \itemize{
    \item Single and multiple predictants
    \item Missing value detection
    \item Data type validation
    \item Sample size requirements
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
    alpha = 0.05,
    Training_data = Training_input,
    X = Predictors,
    Y = Predictants,
    Nmin = 5,
    resolution = 50
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