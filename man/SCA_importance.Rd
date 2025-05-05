\name{SCA_importance}
\alias{SCA_importance}

\title{
	Calculate Variable Importance for a Single SCA Tree
}

\description{
	This function calculates the importance of independent variables in explaining the variability of dependent variables for a single Stepwise Cluster Analysis (SCA) tree using the Wilks' Lambda statistic. The importance is calculated based on the contribution of each variable to the reduction in Wilks' Lambda at each split in the tree.

	For calculating importance scores across all trees in an SCE ensemble, use \code{\link{Wilks_importance}} instead.
}

\usage{
SCA_importance(model)
}

\arguments{
	\item{model}{
		A single SCA tree object containing:
		\itemize{
			\item Tree: Tree structure with Wilks' Lambda values and split information
			\item XName: Names of predictors used
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
		\item Extract Wilks' Lambda values and split information from the tree
		\item Replace negative Wilks' Lambda values with zero
		\item Calculate raw importance for each split:
			\itemize{
				\item Importance = (left_samples + right_samples) / total_samples * (1 - Wilks' Lambda)
			}
		\item Aggregate importance scores by predictor
		\item Normalize importance scores to sum to 1
	}
	
	The function handles:
	\itemize{
		\item Different sets of predictors in the tree
		\item Missing or invalid splits
		\item Both single and multiple predictants
		\item Trees with no splits (returns NULL)
	}

	\strong{Relationship with Wilks_importance:}
	\itemize{
		\item \code{SCA_importance} calculates importance scores for a single SCA tree
		\item \code{\link{Wilks_importance}} calculates importance scores across all trees in an SCE ensemble
		\item Both functions use the same underlying importance calculation method
		\item \code{Wilks_importance} with OOB_weight=FALSE is equivalent to taking the median of \code{SCA_importance} scores across all trees
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
data("Streamflow_training_10var")
data("Streamflow_testing_10var")

## Define the independent (x) and dependent (y) variables
Predictors <- c("Prcp", "SRad", "Tmax", "Tmin", "VP", "smlt", "swvl1", "swvl2", "swvl3", "swvl4")
Predictants <- c("Flow")

## Build a single SCA tree
SCA_tree <- SCA(
	Training_data = Streamflow_training_10var,
	X = Predictors,
	Y = Predictants,
	Nmin = 5,
	alpha = 0.05,
	resolution = 1000
)

## Calculate variable importance for the single tree
Tree_importance <- SCA_importance(SCA_tree)

## Print the results
print("Single tree importance scores:")
print(Tree_importance)

## Visualize the importance scores
barplot(Tree_importance$Importance, 
        names.arg = Tree_importance$col_index,
        main = "Variable Importance (Single Tree)",
        xlab = "Predictors",
        ylab = "Importance Score",
        las = 2)
} 