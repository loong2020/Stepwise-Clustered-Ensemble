\name{Wilks_importance}
\alias{Wilks_importance}
\alias{Wilks_importance}

\title{
	Function to Rank the Importance of Independent Variables
}

\description{
	This function utilizes the Wilks feature importance method to rank the importance of independent variables in explaining the variability in dependent variables. The importance is ranked based on their contribution to the variability of dependent variables.
}

\usage{
Wilks_importance(Model, OOB_weight = FALSE)
}

\references{
Li, Kailong, Guohe Huang, and Brian Baetz. "Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling." Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.
}

\arguments{
	\item{Model}{
		An object representing the trained SCE model, obtained after using the SCE training function.
	}
	\item{OOB_weight}{
		A logical value indicating whether to use out-of-bag (OOB) errors as weights when ranking variable importance. The default value is FALSE, implying that OOB errors are not used in the ranking.
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

## Load the training and testing data files
data("Training_input")
data("Testing_input")

## Define the independent (x) and dependent (y) variables
Predictors <- c("Prcp", "SRad", "Tmax", "Tmin", "VP", "smlt", "swvl1", "swvl2")
Predictants <- c("swvl3", "swvl4")

## Build the SCE model
Model <- SCE(Training_data = Training_input, X = Predictors, Y = Predictants, mfeature = round(0.5 * length(Predictors)), Nmin = 5, Ntree = 48, alpha = 0.05)

## Perform model simulations
Results <- Model_simulation(Testing_data = Testing_input, X = Predictors, model = Model)

## Evaluate the model for specified predictants
Evaluation_swvl3 <- Model_evaluation(Predictant = "swvl3", obs_training = Training_input, obs_testing = Testing_input, sim = Results, Num_predictor = length(Predictors), digits = 2)
Evaluation_swvl4 <- Model_evaluation(Predictant = "swvl4", obs_training = Training_input, obs_testing = Testing_input, sim = Results, Num_predictor = length(Predictors), digits = 2)

## Use the Wilks_importance function to rank the importance of predictors
Importance_score <- Wilks_importance(Model, OOB_weight = FALSE)
}

