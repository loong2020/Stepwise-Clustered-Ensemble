\name{SCE}
\alias{SCE}
\title{
  Stepwise Clustered Ensemble (SCE)
}
\description{
  Builds a Stepwise Clustered Ensemble (SCE) model, which is an ensemble of SCA trees built using bootstrap samples and random feature selection, providing improved prediction accuracy and robustness.
}
\usage{
SCE(Training_data, X, Y, mfeature, Nmin, Ntree, alpha = 0.05, 
    resolution = 1000, verbose = FALSE, parallel = TRUE)
}
\arguments{
  \item{Training_data}{A data.frame containing the training data}
  \item{X}{Character vector of predictor variable names}
  \item{Y}{Character vector of predictant variable names}
  \item{mfeature}{Number of features to randomly select for each tree}
  \item{Nmin}{Minimum number of samples in a leaf node}
  \item{Ntree}{Number of trees in the ensemble}
  \item{alpha}{Significance level for clustering (default: 0.05)}
  \item{resolution}{Resolution for splitting (default: 1000)}
  \item{verbose}{Print progress information (default: FALSE)}
  \item{parallel}{Use parallel processing (default: TRUE)}
}
\value{
  An S3 object of class "SCE" containing the ensemble model.
}
\examples{
\donttest{
  # Load example data
  data(Streamflow_training_10var)
  data(Streamflow_testing_10var)
  
  # Define variables
  Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
  Predictants <- c("Flow")
  
  # Build SCE model
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
  
  # Use S3 methods
  print(sce_model)
  summary(sce_model)
  sce_predictions <- predict(sce_model, Streamflow_testing_10var)
  sce_importance <- importance(sce_model)
  sce_evaluation <- evaluate(sce_model, Streamflow_testing_10var, Streamflow_training_10var)
}
}
\seealso{
  \code{\link{SCA}}, \code{\link{predict}}, \code{\link{importance}}, \code{\link{evaluate}}
} 