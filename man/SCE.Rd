\name{SCE}
\alias{SCE}
\alias{SCE}

\title{
	Build an SCE model
}

\description{
	This function aims to build an SCE model that contains information for prediction and inference
}

\usage{
SCE(Training_data, X, Y, mfeature,Nmin,Ntree,alpha)
}

\arguments{
	\item{Training_data}{
		a data.frame used for training the model.
	}
	\item{X}{
		a string to specify names for independent (x) variable (e.g., c("Prcp","SRad",...,"swvl2")). The length for independent variables must be greater than it for depedent variables.
	}
	\item{Y}{
		a string to specify the name(s) for dependent (y) variable(s) (e.g., c("swvl3","swvl4")).
	}
	\item{mfeature}{
		a value to specify how many features will be drawed for each tree. Default value is round(0.5*length(X)).
	}
	\item{Nmin}{
		a value to specify the minimal number of sample in a leaf node for cutting. Default value is 5.
	}
	\item{Ntree}{
		a value to specify how many trees (snsemble members) will be built. Default value is 48.
	}
	\item{alpha}{
		significance level for clustering, usually in 0.001 - 0.1, default value is 0.05.
	}
}

\author{
Kailong li  <lkl98509509@gmail.com>
}

\references{
Li, Kailong, Guohe Huang, and Brian Baetz. Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling. Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.

Wang, X., G. Huang, Q. Lin, X. Nie, G. Cheng, Y. Fan, Z. Li, Y. Yao, and M. Suo (2013), A stepwise cluster analysis approach for downscaled climate projection - A Canadian case study. Environmental Modelling & Software, 49, 141-151.

Huang, G. (1992). A stepwise cluster analysis method for predicting air quality in an urban environment. Atmospheric Environment (Part B. Urban Atmosphere), 26(3): 349-357.

Liu, Y. Y. and Y. L. Wang (1979). Application of stepwise cluster analysis in medical research. Scientia Sinica, 22(9): 1082-1094.
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
