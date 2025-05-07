\name{SCA_tree_predict}
\alias{SCA_tree_predict}

\title{
	Make Predictions Using a Single SCA Tree
}

\description{
	This function makes predictions using a single Stepwise Cluster Analysis (SCA) tree. It traverses the tree structure based on the predictor values in the test data and returns the predicted values for the predictants. The function includes comprehensive input validation for data types, missing values, and predictor matching.

	For making predictions using an entire SCE ensemble, use \code{\link{Model_simulation}} instead.
}

\usage{
SCA_tree_predict(Testing_data, model)
}

\arguments{
	\item{Testing_data}{
		A data.frame or matrix containing the test data.
		Must include all predictors used in model training.
		Must not contain missing values.
	}
	\item{model}{
		A trained SCA model object returned by the SCA function.
	}
}

\value{
	A list containing:
	\itemize{
		\item A data.frame containing predictions for the test data
	}
}

\details{
	The prediction process involves the following steps:
	\enumerate{
		\item Input validation:
			\itemize{
				\item Data type and structure checks (data.frame or matrix)
				\item Missing value checks
				\item Predictor matching with training data
				\item Numeric data validation
			}
		\item Data preparation:
			\itemize{
				\item Conversion to matrix format
				\item Initialization of prediction matrix
			}
		\item Tree traversal and prediction:
			\itemize{
				\item Processing each test sample through the tree
				\item Generating predictions using leaf node mappings
			}
	}
}

\author{
Kailong Li <lkl98509509@gmail.com>
}

\references{
Li, Kailong, Guohe Huang, and Brian Baetz. "Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling." Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.
}

\examples{
## Load SCE package
library(SCE)
library(parallel)

## Load training and testing data
data("Streamflow_training_10var")
data("Streamflow_testing_10var")

## Define independent (x) and dependent (y) variables
Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
Predictants <- c("Flow")

## Build the SCA model
Model <- SCA(
	Training_data = Streamflow_training_10var,
	X = Predictors,
	Y = Predictants,
	Nmin = 5,
	alpha = 0.05,
	resolution = 1000
)

## Make predictions
Predictions <- SCA_tree_predict(
	Testing_data = Streamflow_testing_10var,
	model = Model
)

} 