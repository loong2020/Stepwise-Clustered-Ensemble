\name{Model_simulation}
\alias{Model_simulation}
\alias{Model_simulation}

\title{
  Perform Simulations using a Trained SCE Model
}

\description{
  This function facilitates the simulation of dependent variables based on the specified independent variables utilizing a previously trained SCE model. It is vital to ensure that the independent variables passed to this function match those used during the model training process.
}

\usage{
Model_simulation(Testing_data, X, model)
}

\arguments{
  \item{Testing_data}{
    A data.frame comprising the data that will be used to test the model. It should contain all the independent variables specified in the X parameter.
  }
  \item{X}{
    A vector of strings specifying the names of independent (x) variables to be used in the simulation (e.g., c("Prcp","SRad",...,"swvl2")). The variable names must exactly match those used during the model training phase.
  }
  \item{model}{
    The trained SCE model object generated through the \code{SCE} function. This model contains the necessary information and parameters for conducting simulation and inferences.
  }
}

\author{
  Kailong Li <lkl98509509@gmail.com>
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

  ## Build the SCE model
  Model <- SCE(Training_data = Training_input, X = Predictors, Y = Predictants, mfeature = round(0.5 * length(Predictors)), Nmin = 5, Ntree = 48, alpha = 0.05)

  ## Conduct simulations using the trained SCE model
  Results <- Model_simulation(Testing_data = Testing_input, X = Predictors, model = Model)
}
