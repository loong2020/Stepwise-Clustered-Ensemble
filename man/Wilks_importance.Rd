\name{Wilks_importance}
\alias{Wilks_importance}

\title{
	Calculate Variable Importance using Wilks' Lambda
}

\description{
	This function calculates the importance of independent variables in explaining the variability of dependent variables using the Wilks' Lambda statistic. The importance is calculated based on the contribution of each variable to the reduction in Wilks' Lambda at each split in the SCA trees. The function supports both unweighted and OOB-weighted importance calculations.
}

\usage{
Wilks_importance(model, OOB_weight = FALSE)
}

\arguments{
	\item{model}{
		A trained SCE model object containing a list of SCA trees. Each tree should contain:
		\itemize{
			\item Tree: Tree structure with Wilks' Lambda values
			\item XName: Names of predictors used
			\item weight: Tree weight (if OOB_weight = TRUE)
		}
	}
	\item{OOB_weight}{
		A logical value indicating whether to weight the importance scores by the tree's OOB performance.
		\itemize{
			\item If TRUE: Importance scores are weighted by each tree's OOB performance
			\item If FALSE: Importance scores are calculated using the median across trees
		}
	}
}

\value{
	A data.frame containing:
	\itemize{
		\item col_index: Names of the predictors
		\item Importance: Normalized importance scores (sum to 1)
	}
}

\details{
	The importance calculation process involves the following steps:
	\enumerate{
		\item Extract Wilks' Lambda values from each tree
		\item Replace negative Wilks' Lambda values with zero
		\item Calculate raw importance for each split:
			\itemize{
				\item Importance = (left_samples + right_samples) / total_samples * (1 - Wilks' Lambda)
			}
		\item Aggregate importance scores by predictor:
			\itemize{
				\item If OOB_weight = TRUE: Weight by tree's OOB performance and sum
				\item If OOB_weight = FALSE: Take median across trees
			}
		\item Normalize importance scores to sum to 1
	}
	
	The function handles:
	\itemize{
		\item Multiple trees in the ensemble
		\item Different sets of predictors in each tree
		\item Missing or invalid splits
		\item Both single and multiple predictants
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
library(hydroGOF)

## Load the training and testing data files
data("Training_input")
data("Testing_input")

## Define the independent (x) and dependent (y) variables
Predictors <- c("Prcp", "SRad", "Tmax", "Tmin", "VP", "smlt", "swvl1", "swvl2")
Predictants <- c("swvl3", "swvl4")

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

## Calculate variable importance without OOB weighting
Importance_unweighted <- Wilks_importance(Model, OOB_weight = FALSE)

## Calculate variable importance with OOB weighting
Importance_weighted <- Wilks_importance(Model, OOB_weight = TRUE)

## Print the results
print("Unweighted importance scores:")
print(Importance_unweighted)

print("OOB-weighted importance scores:")
print(Importance_weighted)
}

