\name{SCE}
\alias{SCE}

\title{
	Build a Stepwise Clustered Ensemble (SCE) Model
}

\description{
	This function builds a Stepwise Clustered Ensemble (SCE) model for multivariate data analysis.
	The SCE model is an ensemble of Stepwise Cluster Analysis (SCA) trees, where each tree is built 
	using bootstrap samples and random feature selection. The function includes comprehensive input 
	validation for data types, missing values, and sample size requirements.
}

\usage{
SCE(Training_data, X, Y, mfeature, Nmin, Ntree, 
    alpha = 0.05, resolution = 1000, verbose = FALSE, parallel = TRUE)
}

\arguments{
	\item{Training_data}{
		A data.frame or matrix containing the training data.
		Must contain all specified predictors and predictants.
		Must not contain missing values.
	}
	\item{X}{
		A character vector specifying the names of independent (predictor) variables
		(e.g., c("Prcp","SRad","Tmax")). Must be present in Training_data.
		All variables must be numeric.
	}
	\item{Y}{
		A character vector specifying the name(s) of dependent (predictant) variable(s)
		(e.g., c("Flow") or c("swvl3","swvl4")). Must be present in Training_data.
		All variables must be numeric.
	}
	\item{mfeature}{
		An integer specifying how many features will be randomly selected for each tree.
		Recommended value is round(0.5 * length(X)).
	}
	\item{Nmin}{
		An integer specifying the minimal number of samples in a leaf node for cutting.
		Must be greater than the number of predictants.
	}
	\item{Ntree}{
		An integer specifying how many trees (ensemble members) will be built.
		Recommended values range from 50 to 500 depending on data complexity.
	}
	\item{alpha}{
		Numeric significance level for clustering, between 0 and 1.
		Default value is 0.05.
	}
	\item{resolution}{
		Numeric value specifying the resolution for splitting.
		Controls the granularity of the search for optimal split points.
		Default value is 1000.
	}
	\item{verbose}{
		A logical value indicating whether to print progress information during model building.
		Default value is FALSE.
	}
	\item{parallel}{
		A logical value indicating whether to use parallel processing for tree construction.
		When TRUE, uses multiple CPU cores for faster computation. When FALSE, processes trees sequentially.
		Default value is TRUE.
	}
}

\value{
	A list containing the ensemble model with the following components:
	\itemize{
		\item \code{Trees}: A list of SCA tree models, each containing:
			\itemize{
				\item \code{Tree}: The SCA tree structure
				\item \code{Map}: Mapping information
				\item \code{XName}: Names of predictors used
				\item \code{YName}: Names of predictants
				\item \code{type}: Mapping type
				\item \code{totalNodes}: Total number of nodes
				\item \code{leafNodes}: Number of leaf nodes
				\item \code{cuttingActions}: Number of cutting actions
				\item \code{mergingActions}: Number of merging actions
				\item \code{OOB_error}: Out-of-bag R-squared error
				\item \code{OOB_sim}: Out-of-bag predictions
				\item \code{Sample}: Bootstrap sample indices
				\item \code{Tree_Info}: Tree-specific information
				\item \code{Training_data}: Training data used for the tree
				\item \code{weight}: Tree weight based on OOB performance
			}
	}
}

\details{
	The SCE model is built using the following steps:
	\enumerate{
		\item \strong{Input Validation}:
			\itemize{
				\item Data type and structure checks
				\item Missing value detection
				\item Numeric data validation
				\item Sample size requirements verification
			}
		\item \strong{Data Preparation}:
			\itemize{
				\item Conversion to appropriate format
				\item Dimension checks
				\item Parameter initialization
			}
		\item \strong{Tree Construction}:
			\itemize{
				\item Generation of bootstrap samples
				\item Random feature selection for each tree
				\item Parallel construction of SCA trees
			}
		\item \strong{Model Evaluation}:
			\itemize{
				\item Calculation of out-of-bag (OOB) errors
				\item Weighting of trees based on OOB performance
			}
	}
	
	The ensemble approach provides improved prediction accuracy and robustness compared to 
	single SCA trees, while the OOB validation provides unbiased performance estimates.
}

\author{
	Kailong Li <lkl98509509@gmail.com>
}

\references{
	Li, Kailong, Guohe Huang, and Brian Baetz. Development of a Wilks feature importance method 
	with improved variable rankings for supporting hydrological inference and modelling. 
	Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.

	Wang, X., G. Huang, Q. Lin, X. Nie, G. Cheng, Y. Fan, Z. Li, Y. Yao, and M. Suo (2013), 
	A stepwise cluster analysis approach for downscaled climate projection - A Canadian case study. 
	Environmental Modelling & Software, 49, 141-151.

	Huang, G. (1992). A stepwise cluster analysis method for predicting air quality in an urban environment. 
	Atmospheric Environment (Part B. Urban Atmosphere), 26(3): 349-357.

	Liu, Y. Y. and Y. L. Wang (1979). Application of stepwise cluster analysis in medical research. 
	Scientia Sinica, 22(9): 1082-1094.
}

\seealso{
	\code{\link{SCA}} for single tree construction,
	\code{\link{SCE_Prediction}} for making predictions,
	\code{\link{Model_simulation}} for comprehensive model evaluation,
	\code{\link{Wilks_importance}} for variable importance analysis,
	\code{\link{RFE_SCE}} for recursive feature elimination
}

\examples{
\donttest{
	## Load required packages
	library(SCE)
	library(parallel)

	## Load example datasets
	data("Streamflow_training_10var")
	data("Streamflow_testing_10var")

	## Define predictors and predictants
	Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
	Predictants <- c("Flow")

	## Build the SCE model
	Model <- SCE(
		Training_data = Streamflow_training_10var,
		X = Predictors,
		Y = Predictants,
		mfeature = round(0.5 * length(Predictors)),
		Nmin = 5,
		Ntree = 48,
		alpha = 0.05,
		resolution = 1000,
		parallel = FALSE
	)

	## Generate predictions for test data
	predictions <- SCE_Prediction(
		X_sample = Streamflow_testing_10var,
		model = Model
	)

	## Conduct comprehensive model evaluation
	Results <- Model_simulation(
		Testing_data = Streamflow_testing_10var,
		model = Model
	)

	## Access different prediction components
	training_predictions <- Results$Training
	validation_predictions <- Results$Validation
	testing_predictions <- Results$Testing

	## Calculate variable importance with OOB weighting (default)
	Importance_weighted <- Wilks_importance(Model)

	## Calculate variable importance without OOB weighting
	Importance_unweighted <- Wilks_importance(Model, OOB_weight = FALSE)

	## Visualize the importance scores
	Importance_ranking_sorted <- Importance_weighted[
		order(-Importance_weighted$Relative_Importance), 
	]
	barplot(
		Importance_ranking_sorted$Relative_Importance,
		names.arg = Importance_ranking_sorted$Predictor,
		las = 2,
		col = "skyblue",
		main = "Variable Importance (SCE)",
		ylab = "Importance",
		xlab = "Predictor"
	)
}
}
