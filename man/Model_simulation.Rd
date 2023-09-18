\name{Model_simulation}
\alias{Model_simulation}
\alias{Model_simulation}

\title{
	Make simulations based on a trained SCE model
}

\description{
	This function aims to simulate dependent variables based on independent ones, given a trained SCE model.
}

\usage{
Model_simulation(Testing_data, X, model)
}

\arguments{
	\item{Testing_data}{
		a data.frame used for testing the model.
	}
	\item{X}{
		a string to specify names for independent (x) variable (e.g., c("Prcp","SRad",...,"swvl2")), must be as same as those used for training.
	}
	\item{model}{
		the trained SCE model.
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

}

