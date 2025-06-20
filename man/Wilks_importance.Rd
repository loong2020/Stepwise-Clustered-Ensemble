\name{Wilks_importance}
\alias{Wilks_importance}

\title{
	Calculate Variable Importance using Wilks' Lambda
}

\description{
	This function calculates the importance of independent variables in explaining the variability of dependent variables using the Wilks' Lambda statistic. The importance is calculated based on the contribution of each variable to the reduction in Wilks' Lambda at each split in the SCA trees. The function supports both unweighted and OOB-weighted importance calculations.

	For calculating importance scores for a single SCA tree, use \code{\link{SCA_importance}} instead.
}

\usage{
Wilks_importance(model, OOB_weight = TRUE)
}

\arguments{
	\item{model}{
		A trained SCE model object containing a list of SCA trees. Each tree should contain:
		\itemize{
			\item Tree: Tree structure with Wilks' Lambda values and split information
			\item XName: Names of predictors used
			\item weight: Tree weight (if OOB_weight = TRUE)
		}
	}
	\item{OOB_weight}{
		A logical value indicating whether to weight the importance scores by the tree's OOB performance.
		\itemize{
			\item If TRUE (default): Importance scores are weighted by each tree's OOB performance
			\item If FALSE: Importance scores are calculated using the median across trees
		}
	}
}

\value{
	A data.frame containing:
	\itemize{
		\item Predictor: Names of the predictors
		\item Relative_Importance: Normalized importance scores (sum to 1)
	}
}

\details{
	The importance calculation process involves the following steps:
	\enumerate{
		\item Extract Wilks' Lambda values and split information from each tree
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
		\item Trees with no splits (returns NULL for those trees)
	}

	\strong{Relationship with SCA_importance:}
	\itemize{
		\item \code{Wilks_importance} calculates importance scores across all trees in an SCE ensemble
		\item \code{\link{SCA_importance}} calculates importance scores for a single SCA tree
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

\seealso{
  \code{\link{SCE}}
}

