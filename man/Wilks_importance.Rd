\name{Wilks_importance}
\alias{Wilks_importance}
\alias{Wilks_importance}

\title{
	A function for importance ranking
}

\description{
	This function aims to rank the importance of independent variables in explaining the variability of dependent variables.
}

\usage{
Wilks_importance(Model,OOB_weight=FALSE)
}

\references{
Li, Kailong, Guohe Huang, and Brian Baetz. Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling. Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.


\arguments{
	\item{Model}{
		the trained SCE model.
	}
	\item{OOB_weight}{
		a logical value to specify wheather to use out-of-bag errors as weights for variable ranking. Default value is FALSE.
	}
}

\author{
Kailong li  <lkl98509509@gmail.com>
}

\examples{
## Load SCE package
library(SCE)

## Training data file
Training_input <- data("Training_input")
## Testing data file
Testing_input <- data("Testing_input")

## Define independent (x) and dependent (y) variables
Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2")
Predictants <- c("swvl3","swvl4")

## Build the SCE model
Model <- SCE(Training_data=Training_input, X=Predictors, Y=Predictants, mfeature=round(0.5*length(Predictors)),Nmin=5,Ntree=48,alpha=0.05)

## Make predictions
Results <- Model_simulation(Testing_data=Testing_input, X=Predictors, model=Model)

## Evaluate the model over a specified predictant
Evaluation_swvl3 <- Model_evaluation(Predictant="swvl3",obs_training=Training_input,obs_testing=Testing_input,sim=Results,Num_predictor=length(Predictors),digits=2)

Evaluation_swvl4 <- Model_evaluation(Predictant="swvl4",obs_training=Training_input,obs_testing=Testing_input,sim=Results,Num_predictor=length(Predictors),digits=2)

## Rank the importance of predictors
Importance_score <- Wilks_importance(Model,OOB_weight=FALSE)
}
