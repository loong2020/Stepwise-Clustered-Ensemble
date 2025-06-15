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
\dontrun{
  ## Load SCE package and the supporting packages
  library(SCE)
  library(parallel)

  data(Streamflow_training_22var)
  data(Streamflow_testing_22var)

  # Define predictors and predictants
  Predictors <- c(
    "Precipitation", "Radiation", "Tmax", "Tmin", "VP",
    "Precipitation_2Mon", "Radiation_2Mon", "Tmax_2Mon", "Tmin_2Mon", "VP_2Mon",
    "PNA", "Nino3.4", "IPO", "PDO",
    "PNA_lag1", "Nino3.4_lag1", "IPO_lag1", "PDO_lag1",
    "PNA_lag2", "Nino3.4_lag2", "IPO_lag2", "PDO_lag2"
  )
  Predictants <- c("Flow")

  # Perform RFE
  set.seed(123)
  result <- RFE_SCE(
    Training_data = Streamflow_training_22var,
    Testing_data = Streamflow_testing_22var,
    Predictors = Predictors,
    Predictant = Predictants,
    Nmin = 5,
    Ntree = 48,
    alpha = 0.05,
    resolution = 1000,
    step = 3  # Number of predictors to remove at each iteration
  )

  ## Access results
  summary <- RFE_results$summary
  performances <- RFE_results$performances
  importance_scores <- RFE_results$importance_scores

  ## Print summary
  print("RFE Summary:")
  print(summary)
}
} 