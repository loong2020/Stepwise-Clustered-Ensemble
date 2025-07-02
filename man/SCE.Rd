\name{SCE}
\alias{SCE}
\alias{SCA}

\title{
	Stepwise Clustered Ensemble (SCE) and Stepwise Cluster Analysis (SCA) Models
}

\description{
	This package provides two main modeling approaches:
	
	\strong{SCA (Stepwise Cluster Analysis):} A single tree model that recursively partitions the data space based on Wilks' Lambda statistic, creating a tree structure for prediction.
	
	\strong{SCE (Stepwise Clustered Ensemble):} An ensemble of SCA trees built using bootstrap samples and random feature selection, providing improved prediction accuracy and robustness.
	
	Both functions include comprehensive input validation for data types, missing values, and sample size requirements, and support both single and multiple predictants.
}

\usage{
SCA(Training_data, X, Y, Nmin, alpha = 0.05, resolution = 1000, verbose = FALSE)

SCE(Training_data, X, Y, mfeature, Nmin, Ntree,
	alpha = 0.05, resolution = 1000, verbose = FALSE, parallel = TRUE)
}

\arguments{
	\item{Training_data}{
		A data.frame or matrix containing the training data.
		Must include all specified predictors and predictants.
		Must not contain missing values.
	}
	\item{X}{
		A character vector specifying the names of independent variables
		(e.g., c("Prcp","SRad","Tmax")). Must be present in Training_data.
		All variables must be numeric.
	}
	\item{Y}{
		A character vector specifying the name(s) of dependent variable(s)
		(e.g., c("swvl3","swvl4")). Must be present in Training_data.
		All variables must be numeric.
	}
	\item{Nmin}{
		Integer specifying the minimal number of samples in a leaf node for cutting.
		Must be a positive number and less than the sample size.
	}
	\item{mfeature}{
		An integer specifying how many features will be randomly selected for each tree.
		Recommended value is round(0.5 * length(X)). Only used for SCE.
	}
	\item{Ntree}{
		An integer specifying how many trees (ensemble members) will be built.
		Recommended values range from 50 to 500 depending on data complexity. Only used for SCE.
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
		Default value is TRUE. Only used for SCE.
	}
}

\value{
	\strong{For SCA:} An S3 object of class "SCA" containing:
	\itemize{
		\item Tree: The SCA tree structure
		\item Map: Mapping information for predictions
		\item XName: Names of predictors used
		\item YName: Names of predictants
		\item type: Mapping type (currently "mean")
		\item totalNodes: Total number of nodes in the tree
		\item leafNodes: Number of leaf nodes
		\item cuttingActions: Number of cutting actions performed
		\item mergingActions: Number of merging actions performed
		\item call: Function call
	}
	
	\strong{For SCE:} An S3 object of class "SCE" containing the ensemble model with the following components:
	\itemize{
		\item \code{trees}: A list of SCA tree models, each containing:
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
		\item \code{predictors}: Names of predictor variables
		\item \code{predictants}: Names of predictant variables
		\item \code{parameters}: Model parameters
		\item \code{call}: Function call
	}
	
	Both objects support S3 methods: \code{print()}, \code{summary()}, 
	\code{predict()}, \code{importance()}, and \code{evaluate()}.
}

\details{
	\strong{Model Building Process:}
	
	\emph{SCA (Single Tree):}
	\enumerate{
		\item Input validation (data types, missing values, sample size requirements)
		\item Data preparation (conversion to matrix format, parameter initialization)
		\item Tree construction (recursive partitioning based on Wilks' Lambda)
	}
	
	\emph{SCE (Ensemble):}
	\enumerate{
		\item Input validation (data types, missing values, sample size requirements)
		\item Data preparation (conversion to matrix format, parameter initialization)
		\item Tree construction (bootstrap samples, random feature selection, parallel SCA tree building)
		\item Model evaluation (OOB error calculation, tree weighting)
	}
	
	\strong{Key Differences:}
	\itemize{
		\item \strong{SCA}: Single tree, deterministic, faster training, potentially less robust
		\item \strong{SCE}: Multiple trees, ensemble approach, improved accuracy, OOB validation, parallel processing
	}
	
	\strong{When to Use:}
	\itemize{
		\item \strong{SCA}: Quick exploration, simple relationships, limited computational resources
		\item \strong{SCE}: Production models, complex relationships, when accuracy is critical
	}
}

\author{
	Xiuquan Wang <xxwang@upei.ca> (original SCA)
	Kailong Li <lkl98509509@gmail.com> (Resolution-search-based SCA and SCE ensemble)
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
	\code{\link{predict}}, \code{\link{importance}}, \code{\link{evaluate}} for S3 methods,
	\code{\link{RFE_SCE}} for recursive feature elimination
}

\examples{
\donttest{
	## Load SCE package
	library(SCE)

	## Load training and testing data
	data("Streamflow_training_10var")
	data("Streamflow_testing_10var")

	## Define independent (x) and dependent (y) variables
	Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
	Predictants <- c("Flow")

	## Example 1: Build SCA model (single tree)
	sca_model <- SCA(
		Training_data = Streamflow_training_10var,
		X = Predictors,
		Y = Predictants,
		Nmin = 5,
		alpha = 0.05,
		resolution = 1000
	)
	
	## Use S3 methods for SCA model inspection
	print(sca_model)
	summary(sca_model)
	
	## Make predictions using S3 method
	sca_predictions <- predict(sca_model, Streamflow_testing_10var)
	
	## Calculate variable importance using S3 method
	sca_importance <- importance(sca_model)
	
	## Evaluate SCA model performance using S3 method
	sca_evaluation <- evaluate(
		object = sca_model,
		Testing_data = Streamflow_testing_10var,
		Predictant = Predictants
	)

	## Example 2: Build SCE model (ensemble)
	sce_model <- SCE(
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

	## Use S3 methods for SCE model inspection
	print(sce_model)
	summary(sce_model)

	## Generate predictions using S3 method
	sce_predictions <- predict(sce_model, Streamflow_testing_10var)

	## Access different prediction components
	training_predictions <- sce_predictions$Training
	validation_predictions <- sce_predictions$Validation
	testing_predictions <- sce_predictions$Testing

	## Calculate variable importance using S3 method
	sce_importance <- importance(sce_model)

	## Evaluate SCE model performance using S3 method
	sce_evaluation <- evaluate(
		object = sce_model,
		Testing_data = Streamflow_testing_10var,
		Training_data = Streamflow_training_10var,
		Predictant = Predictants
	)
}
}
