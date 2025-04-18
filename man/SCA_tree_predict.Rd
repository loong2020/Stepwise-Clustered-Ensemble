\name{SCA_tree_predict}
\alias{SCA_tree_predict}

\title{
	Make Predictions Using a Single SCA Tree
}

\description{
	This function makes predictions using a single Stepwise Cluster Analysis (SCA) tree. It traverses the tree structure based on the predictor values in the test data and returns the predicted values for the predictants.

	For making predictions using an entire SCE ensemble, use \code{\link{Model_simulation}} instead.
}

\usage{
SCA_tree_predict(Test_data, X, model)
}

\arguments{
	\item{Test_data}{
		A data.frame containing the test data with the same predictor variables as used in training.
	}
	\item{X}{
		A character vector containing the names of the predictor variables.
	}
	\item{model}{
		A single SCA tree object containing:
		\itemize{
			\item Tree: Tree structure with split information
			\item XName: Names of predictors used
			\item YName: Names of predictants
		}
	}
}

\value{
	A matrix containing the predicted values for each predictant. Each column corresponds to a predictant, and each row corresponds to a sample in the test data.
}

\details{
	The prediction process involves the following steps:
	\enumerate{
		\item For each sample in the test data:
			\itemize{
				\item Start at the root node of the tree
				\item Traverse down the tree based on the split conditions
				\item When reaching a leaf node, use its predictant values as the prediction
			}
		\item Return all predictions as a matrix
	}
	
	The function handles:
	\itemize{
		\item Both single and multiple predictants
		\item Missing values in the test data
		\item Different sets of predictors in the tree
		\item Invalid or missing splits
	}

	\strong{Relationship with Model_simulation:}
	\itemize{
		\item \code{SCA_tree_predict} makes predictions using a single SCA tree
		\item \code{\link{Model_simulation}} makes predictions using the entire SCE ensemble
		\item \code{Model_simulation} internally uses \code{SCA_tree_predict} for each tree in the ensemble
	}
}

\author{
Kailong Li <lkl98509509@gmail.com>
}

\references{
Li, Kailong, Guohe Huang, and Brian Baetz. "Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling." Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.
}

\examples{
## Load SCE package and the supporting packages
library(SCE)
library(parallel)

## Load the training and testing data files
data("Training_input")
data("Testing_input")

## Define the independent (x) and dependent (y) variables
Predictors <- c("Prcp", "SRad", "Tmax", "Tmin", "VP", "smlt", "swvl1", "swvl2")
Predictants <- c("swvl3", "swvl4")

## Build a single SCA tree
SCA_tree <- SCA(
	Training_data = Training_input,
	X = Predictors,
	Y = Predictants,
	Nmin = 5,
	alpha = 0.05,
	resolution = 1000
)

## Make predictions using the single tree
Predictions <- SCA_tree_predict(
	Test_data = Testing_input,
	X = Predictors,
	model = SCA_tree
)

## Print the first few predictions
print("First few predictions:")
print(head(Predictions))

## Evaluate the predictions
Evaluation <- SCA_Model_evaluation(
	Testing_data = Testing_input,
	Simulations = Predictions,
	Predictant = Predictants,
	digits = 2
)

print("Model evaluation results:")
print(Evaluation)
} 