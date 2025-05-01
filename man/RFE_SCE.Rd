\name{RFE_SCE}
\alias{RFE_SCE}

\title{
  Recursive Feature Elimination for SCE Models
}

\description{
  This function implements Recursive Feature Elimination (RFE) to identify the most important predictors for SCE models. It iteratively removes the least important predictors based on Wilks' feature importance scores and evaluates model performance. The function supports both single and multiple predictants, with comprehensive input validation and performance tracking across iterations.
}

\usage{
RFE_SCE(
  Training_data,
  Testing_data,
  Predictors,
  Predictant,
  Nmin,
  Ntree,
  mfeature,
  alpha = 0.05,
  resolution = 1000,
  step = 1
)
}

\arguments{
  \item{Training_data}{
    A data.frame containing the training data. Must include all specified predictors and predictants.
  }
  \item{Testing_data}{
    A data.frame containing the testing data. Must include all specified predictors and predictants.
  }
  \item{Predictors}{
    A character vector specifying the names of independent variables to be evaluated
    (e.g., c("Prcp","SRad","Tmax")). Must contain at least 2 elements.
  }
  \item{Predictant}{
    A character vector specifying the name(s) of dependent variable(s)
    (e.g., c("swvl3","swvl4")). Must be non-empty.
  }
  \item{Nmin}{
    Integer specifying the minimal number of samples in a leaf node for cutting.
  }
  \item{Ntree}{
    Integer specifying the number of trees in the ensemble.
  }
  \item{mfeature}{
    Integer specifying how many features will be randomly selected for each tree.
  }
  \item{alpha}{
    Numeric significance level for clustering, between 0 and 1.
    Default value is 0.05.
  }
  \item{resolution}{
    Numeric value specifying the resolution for splitting.
    Default value is 1000.
  }
  \item{step}{
    Integer specifying the number of predictors to remove at each iteration.
    Must be between 1 and (number of predictors - number of predictants).
    Default value is 1.
  }
}

\value{
  A list containing:
  \itemize{
    \item summary: Data.frame with columns:
      \itemize{
        \item n_predictors: Number of predictors at each iteration
        \item predictors: Comma-separated list of predictors used
      }
    \item performances: List of performance evaluations for each iteration
      \itemize{
        \item For single predictant: Direct performance data.frame
        \item For multiple predictants: Named list of performance data.frames
      }
    \item importance_scores: List of Wilks' importance scores for each iteration
  }
}

\details{
  The RFE process involves the following steps:
  \enumerate{
    \item Input validation:
      \itemize{
        \item Data frame structure validation
        \item Predictor and predictant validation
        \item Step size validation
      }
    \item Initialization:
      \itemize{
        \item Set up history tracking structures
        \item Initialize current predictor set
      }
    \item Main RFE loop (continues while predictors > predictants + 2):
      \itemize{
        \item Train SCE model with current predictors
        \item Generate predictions using Model_simulation
        \item Evaluate model using SCE_Model_evaluation
        \item Store performance metrics and importance scores
        \item Remove least important predictors based on Wilks' scores
      }
  }
  
  The function handles:
  \itemize{
    \item Single and multiple predictants
    \item Performance tracking across iterations
    \item Importance score calculation
    \item Step-wise predictor removal
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

  ## Perform RFE
  RFE_results <- RFE_SCE(
    Training_data = Training_input,
    Testing_data = Testing_input,
    Predictors = Predictors,
    Predictant = Predictants,
    Nmin = 5,
    Ntree = 40,
    mfeature = round(0.5 * length(Predictors)),
    alpha = 0.05,
    resolution = 1000,
    step = 1
  )

  ## Access results
  summary <- RFE_results$summary
  performances <- RFE_results$performances
  importance_scores <- RFE_results$importance_scores

  ## Print summary
  print("RFE Summary:")
  print(summary)
} 