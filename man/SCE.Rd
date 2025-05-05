\name{SCE}
\alias{SCE}

\title{
	Build a Stepwise Clustered Ensemble (SCE) model
}

\description{
	This function builds an SCE model that contains information for prediction and inference. The SCE model is an ensemble of Stepwise Cluster Analysis (SCA) trees, where each tree is built using bootstrap samples and random feature selection. The function includes comprehensive input validation for data types, missing values, and sample size requirements.
}

\usage{
SCE(Training_data, X, Y, mfeature, Nmin, Ntree, alpha = 0.05, resolution = 1000)
}

\arguments{
	\item{Training_data}{
		A data.frame or matrix used for training the model.
		Must contain all specified predictors and predictants.
		Must not contain missing values.
	}
	\item{X}{
		A character vector specifying names for independent (predictor) variables
		(e.g., c("Prcp","SRad","Tmax")). Must be present in Training_data.
		All variables must be numeric.
	}
	\item{Y}{
		A character vector specifying the name(s) for dependent (predictant) variable(s)
		(e.g., c("swvl3","swvl4")). Must be present in Training_data.
		All variables must be numeric.
	}
	\item{mfeature}{
		An integer specifying how many features will be randomly selected for each tree.
		Recommended value is round(0.5*length(X)).
	}
	\item{Nmin}{
		An integer specifying the minimal number of samples in a leaf node for cutting.
		Must be greater than the number of predictants.
	}
	\item{Ntree}{
		An integer specifying how many trees (ensemble members) will be built.
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
}

\value{
	A list containing:
	\itemize{
		\item A list of SCA tree models, each containing:
			\itemize{
				\item Tree: The SCA tree structure
				\item Map: Mapping information
				\item XName: Names of predictors used
				\item YName: Names of predictants
				\item type: Mapping type
				\item totalNodes: Total number of nodes
				\item leafNodes: Number of leaf nodes
				\item cuttingActions: Number of cutting actions
				\item mergingActions: Number of merging actions
				\item OOB_error: Out-of-bag R-squared error
				\item OOB_sim: Out-of-bag predictions
				\item Sample: Bootstrap sample indices
				\item Tree_Info: Tree-specific information
				\item Training_data: Training data used for the tree
				\item weight: Tree weight based on OOB performance
			}
	}
}

\details{
	The SCE model is built using the following steps:
	\enumerate{
		\item Input validation:
			\itemize{
				\item Data type and structure checks
				\item Missing value checks
				\item Numeric data validation
				\item Sample size requirements
			}
		\item Data preparation:
			\itemize{
				\item Conversion to matrix format
				\item Dimension checks
				\item Parameter initialization
			}
		\item Tree construction:
			\itemize{
				\item Generation of bootstrap samples
				\item Random feature selection
				\item Parallel construction of SCA trees
			}
		\item Model evaluation:
			\itemize{
				\item Calculation of out-of-bag (OOB) errors
				\item Weighting of trees based on OOB performance
			}
	}
}

\author{
	Kailong Li <lkl98509509@gmail.com>
}

\references{
	Li, Kailong, Guohe Huang, and Brian Baetz. Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling. Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.

	Wang, X., G. Huang, Q. Lin, X. Nie, G. Cheng, Y. Fan, Z. Li, Y. Yao, and M. Suo (2013), A stepwise cluster analysis approach for downscaled climate projection - A Canadian case study. Environmental Modelling & Software, 49, 141-151.

	Huang, G. (1992). A stepwise cluster analysis method for predicting air quality in an urban environment. Atmospheric Environment (Part B. Urban Atmosphere), 26(3): 349-357.

	Liu, Y. Y. and Y. L. Wang (1979). Application of stepwise cluster analysis in medical research. Scientia Sinica, 22(9): 1082-1094.
}

\examples{
	## Load SCE package and the supporting packages
	library(SCE)
	library(parallel)

	## Training data file
	data("Streamflow_training_10var")
	## Testing data file
	data("Streamflow_testing_10var")

	## Define independent (x) and dependent (y) variables
	Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
	Predictants <- c("Flow")

	## Build the SCE model
	Model <- SCE(
		Training_data = Streamflow_training_10var,
		X = Predictors,
		Y = Predictants,
		mfeature = round(0.5*length(Predictors)),
		Nmin = 5,
		Ntree = 48,
		alpha = 0.05,
		resolution = 1000
	)
}
