\name{Model_evaluation}
\alias{Model_evaluation}
\alias{Model_evaluation}

\title{
	Evaluate SCE model performances
}

\description{
	This function aims to evaluate SCE model performances over training, out-of-bag validation, and testing datasets. The function outputs 22 evaluation metrics.
}

\usage{
Model_evaluation(Predictant,obs_training,obs_testing,sim,Num_predictor,digits)
}

\arguments{
	\item{Predictant}{
		a string to specify the name for ONE of dependent (y) variables (e.g., "streamflow").
	}
	\item{obs_training}{
		the data.frame used for training the model.
	}
	\item{obs_testing}{
		the data.frame used for testing the model.
	}
	\item{sim}{
		a list including model simulations for training, out-of-bag validation, and testing.
	}
	\item{Num_predictor}{
		the number of independent (x) variables used for calculating the ajusted R-squared.
	}
	\item{digits}{
		the number of digits to keep for evaluation metrics.
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

}

