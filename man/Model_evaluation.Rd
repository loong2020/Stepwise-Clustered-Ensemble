\name{Model_evaluation}
\alias{Model_evaluation}
\alias{SCE_Model_evaluation}
\alias{SCA_Model_evaluation}

\title{
  Evaluate SCE and SCA Model Performances
}

\description{
  These functions facilitate the evaluation of SCE and SCA model performances over training,
  out-of-bag validation, and testing datasets, providing a comprehensive summary
  through 6 distinct evaluation metrics: MAE, RMSE, NSE, Log-NSE, Adjusted R-squared, and KGE.
}

\usage{
SCE_Model_evaluation(
  Testing_data,
  Training_data,
  Simulations,
  Predictant,
  Num_predictor,
  digits = 3
)

SCA_Model_evaluation(
  Testing_data,
  Simulations,
  Predictant,
  Num_predictor,
  digits = 3
)
}

\arguments{
  \item{Testing_data}{
    A data.frame comprising the observations used during the model testing phase.
    Must contain all specified predictants.
  }
  \item{Training_data}{
    A data.frame comprising the observations used during the model training phase.
    Required only for SCE_Model_evaluation.
  }
  \item{Simulations}{
    A list containing model simulations:
    \itemize{
      \item For SCE: must contain 'Training', 'Validation', and 'Testing' components
      \item For SCA: must contain 'Testing_sim' component
    }
    The structure should align with the output generated by the respective model training function.
  }
  \item{Predictant}{
    A character vector specifying the name(s) of the dependent (y) variable(s) to be evaluated
    (e.g., c("swvl3", "swvl4")). The specified names must exactly match those used in model training.
  }
  \item{Num_predictor}{
    An integer specifying the number of independent (x) variables utilized in the model.
    This value is used to calculate adjusted R-squared metrics and should match the number
    of independent variables specified during model training.
  }
  \item{digits}{
    An integer specifying the number of decimal places to retain when reporting
    evaluation metrics. Default value is 3.
  }
}

\value{
  A data.frame or list of data.frames containing the following metrics:
  \itemize{
    \item MAE: Mean Absolute Error
    \item RMSE: Root Mean Square Error
    \item NSE: Nash-Sutcliffe Efficiency
    \item Log.NSE: Log-transformed Nash-Sutcliffe Efficiency
    \item Adj_R2: Adjusted R-squared
    \item kge: Kling-Gupta Efficiency
  }
  
  For SCE_Model_evaluation:
  \itemize{
    \item If single predictant: Returns a data.frame with columns "Training", "Validation", and "Testing"
    \item If multiple predictants: Returns a list of data.frames, one for each predictant
  }
  
  For SCA_Model_evaluation:
  \itemize{
    \item If single predictant: Returns a data.frame with column "Testing"
    \item If multiple predictants: Returns a list of data.frames, one for each predictant
  }
}

\details{
  The evaluation metrics are calculated as follows:
  \enumerate{
    \item MAE: mean(abs(obs - sim))
    \item RMSE: sqrt(mean((obs - sim)^2))
    \item NSE: 1 - (sum((obs - sim)^2) / sum((obs - mean(obs))^2))
    \item Log.NSE: NSE calculated on log-transformed values
    \item Adj_R2: 1 - (1 - R2) * (n-1)/(n-p-1), where n is sample size and p is number of predictors
    \item KGE: 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2), where r is correlation, alpha is ratio of standard deviations, and beta is ratio of means
  }
  
  Input validation includes checks for:
  \enumerate{
    \item Data frame structure of input data
    \item Presence of required components in Simulations list
    \item Existence of predictants in both data and simulations
    \item Matching row counts between data and simulations
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

  ## Training data file
  data("Training_input")
  ## Testing data file
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

  ## Make predictions
  Results <- Model_simulation(
    Testing_data = Testing_input,
    Training_data = Training_input,
    model = Model
  )

  ## Evaluate the SCE model
  Evaluation <- SCE_Model_evaluation(
    Testing_data = Testing_input,
    Training_data = Training_input,
    Simulations = Results,
    Predictant = Predictants,
    Num_predictor = length(Predictors),
    digits = 3
  )

  ## Evaluate the SCA model
  SCA_Evaluation <- SCA_Model_evaluation(
    Testing_data = Testing_input,
    Simulations = Results,
    Predictant = Predictants,
    Num_predictor = length(Predictors),
    digits = 3
  )
}
