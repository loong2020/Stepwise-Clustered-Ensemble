\name{RFE_SCE}
\alias{RFE_SCE}

\title{
  Recursive Feature Elimination for SCE Models
}

\description{
  This function implements Recursive Feature Elimination (RFE) to identify the most important predictors for SCE models. It iteratively removes the least important predictors based on Wilks' feature importance scores and evaluates model performance using various metrics. The function supports both single and multiple predictants.
}

\usage{
RFE_SCE(
  Training_data,
  Testing_data,
  Predictors,
  Predictant,
  alpha = 0.05,
  Nmin = 5,
  Ntree = 40,
  resolution = 50,
  metric = "rmse",
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
  \item{alpha}{
    Numeric significance level for clustering, usually between 0.001 and 0.1.
    Default value is 0.05.
  }
  \item{Nmin}{
    Integer specifying the minimal number of samples in a leaf node for cutting.
    Default value is 5.
  }
  \item{Ntree}{
    Integer specifying the number of trees in the ensemble.
    Default value is 40.
  }
  \item{resolution}{
    Numeric value specifying the resolution for splitting.
    Default value is 50.
  }
  \item{metric}{
    Character string specifying the performance metric to optimize:
    \itemize{
      \item For single predictant: "rmse", "mae", "nse", "log_nse", "Adj_R2", "kge"
      \item For multiple predictants: "nse", "log_nse", "Adj_R2", "kge"
    }
    Default value is "rmse".
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
    \item best_predictors: Character vector of the optimal set of predictors
    \item best_metric: Numeric value of the best performance metric achieved
    \item history: List containing:
      \itemize{
        \item summary: Data.frame with columns:
          \itemize{
            \item n_predictors: Number of predictors at each iteration
            \item predictors: Comma-separated list of predictors used
            \item current_metric: Performance metric value
          }
        \item performances: List of performance evaluations for each iteration
      }
    \item final_model: The SCE model trained with the best set of predictors
  }
}

\details{
  The RFE process involves the following steps:
  \enumerate{
    \item Input validation and initialization
    \item Iterative process:
      \itemize{
        \item Train SCE model with current predictors
        \item Make predictions on testing data
        \item Evaluate model performance
        \item Calculate predictor importance using Wilks' method
        \item Remove least important predictors
      }
    \item Track performance history
    \item Return best performing predictor set
  }
  
  The function handles:
  \itemize{
    \item Single and multiple predictants
    \item Different performance metrics
    \item Step-wise predictor removal
    \item Performance tracking across iterations
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
    alpha = 0.05,
    Nmin = 5,
    Ntree = 40,
    resolution = 50,
    metric = "nse",
    step = 1
  )

  ## Access results
  best_predictors <- RFE_results$best_predictors
  best_metric <- RFE_results$best_metric
  performance_history <- RFE_results$history
  final_model <- RFE_results$final_model

  ## Print summary
  print("Best predictors:")
  print(best_predictors)
  print("Best metric value:")
  print(best_metric)
} 