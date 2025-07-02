\name{RFE_SCE}
\alias{RFE_SCE}
\alias{Plot_RFE}

\title{
  Recursive Feature Elimination for SCE Models
}

\description{
  This function implements Recursive Feature Elimination (RFE) to identify the most important predictors for SCE models. It iteratively removes the least important predictors based on Wilks' feature importance scores and evaluates model performance. The function supports both single and multiple predictants, with comprehensive input validation and performance tracking across iterations.
  
  The package also provides a \code{Plot_RFE} function for visualizing RFE results, showing validation and testing R² values as a function of the number of predictors.
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
  step = 1,
  verbose = TRUE,
  parallel = TRUE
)

Plot_RFE(
  rfe_result,
  main = "Validation and Testing R² vs Number of Predictors",
  col_validation = "blue",
  col_testing = "red",
  pch = 16,
  lwd = 2,
  cex = 1.2,
  legend_pos = "bottomleft",
  ...
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
  \item{verbose}{
    A logical value indicating whether to print progress information during RFE iterations.
    Default value is TRUE.
  }
  \item{parallel}{
    A logical value indicating whether to use parallel processing for SCE model construction.
    When TRUE, uses multiple CPU cores for faster computation. When FALSE, processes trees sequentially.
    Default value is TRUE.
  }
  
  \strong{Plot_RFE Arguments:}
  \item{rfe_result}{
    The result object from RFE_SCE function containing summary and performances components.
  }
  \item{main}{
    Title for the plot. Default is "Validation and Testing R² vs Number of Predictors".
  }
  \item{col_validation}{
    Color for validation line. Default is "blue".
  }
  \item{col_testing}{
    Color for testing line. Default is "red".
  }
  \item{pch}{
    Point character for markers. Default is 16 (filled circle).
  }
  \item{lwd}{
    Line width. Default is 2.
  }
  \item{cex}{
    Point size. Default is 1.2.
  }
  \item{legend_pos}{
    Position of legend. Default is "bottomleft".
  }
  \item{...}{
    Additional arguments passed to plot function.
  }
}

\value{
  \strong{RFE_SCE:} A list containing:
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
  
  \strong{Plot_RFE:} Invisibly returns a list containing:
  \itemize{
    \item n_predictors: Vector of predictor counts
    \item validation_r2: Vector of validation R² values
    \item testing_r2: Vector of testing R² values
  }
}

\details{
  \strong{RFE_SCE Process:}
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
  
  \strong{Plot_RFE Function:}
  Creates a base R plot showing validation and testing R² values as a function of the number of predictors during the RFE process. The function:
  \itemize{
    \item Extracts R² values from RFE results
    \item Converts formatted strings to numeric values
    \item Creates a line plot with points and lines
    \item Includes a legend distinguishing validation and testing performance
    \item Supports customization of colors, line styles, and plot appearance
    \item Uses only base R graphics (no external dependencies)
  }
}

\author{
  Kailong Li <lkl98509509@gmail.com>
}

\examples{
\donttest{
#   # This example is computationally intensive and may take a long time to run.
#   # It is recommended to run this example on a machine with a high-performance CPU.
# 
#   ## Load SCE package and the supporting packages
#   library(SCE)
#   library(parallel)
# 
#   data(Streamflow_training_22var)
#   data(Streamflow_testing_22var)
# 
#   # Define predictors and predictants
#   Predictors <- c(
#     "Precipitation", "Radiation", "Tmax", "Tmin", "VP",
#     "Precipitation_2Mon", "Radiation_2Mon", "Tmax_2Mon", "Tmin_2Mon", "VP_2Mon",
#     "PNA", "Nino3.4", "IPO", "PDO",
#     "PNA_lag1", "Nino3.4_lag1", "IPO_lag1", "PDO_lag1",
#     "PNA_lag2", "Nino3.4_lag2", "IPO_lag2", "PDO_lag2"
#   )
#   Predictants <- c("Flow")
# 
#   # Perform RFE
#   set.seed(123)
#   result <- RFE_SCE(
#     Training_data = Streamflow_training_22var,
#     Testing_data = Streamflow_testing_22var,
#     Predictors = Predictors,
#     Predictant = Predictants,
#     Nmin = 5,
#     Ntree = 48,
#     alpha = 0.05,
#     resolution = 1000,
#     step = 3,  # Number of predictors to remove at each iteration
#     verbose = TRUE,
#     parallel = TRUE
#   )
#
#   ## Access results
#   summary <- result$summary
#   performances <- result$performances
#   importance_scores <- result$importance_scores
#
#   ## Plot RFE results
#   Plot_RFE(result)
#
#   ## Customized plot
#   Plot_RFE(result, 
#            main = "My RFE Results",
#            col_validation = "darkblue",
#            col_testing = "darkred",
#            lwd = 3,
#            cex = 1.5)
#
#   ## Note: The RFE_SCE function internally uses S3 methods for SCE models
#   ## including importance() and evaluate() for model analysis
# 
# }
} 

\seealso{
  See the generic functions \code{importance} and \code{evaluate} for SCE objects.
  For visualization of RFE results, see \code{Plot_RFE}.
} 