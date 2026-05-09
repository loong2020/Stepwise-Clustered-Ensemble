\name{sce}
\alias{sce}
\alias{SCE}
\title{
  Stepwise Clustered Ensemble (SCE)
}
\description{
  Builds a Stepwise Clustered Ensemble (SCE) model, which is an ensemble of SCA trees built using bootstrap samples and random feature selection, providing improved prediction accuracy and robustness.
}
\usage{
sce(training_data, x, y, mfeature, nmin, ntree, alpha = 0.05,
    resolution = 1000, verbose = FALSE, parallel = TRUE)
}
\arguments{
  \item{training_data}{A data.frame containing the training data}
  \item{x}{Character vector of predictor variable names}
  \item{y}{Character vector of predictant variable names}
  \item{mfeature}{Number of features to randomly select for each tree}
  \item{nmin}{Minimum number of samples in a leaf node}
  \item{ntree}{Number of trees in the ensemble}
  \item{alpha}{Significance level for clustering (default: 0.05)}
  \item{resolution}{Resolution for splitting (default: 1000)}
  \item{verbose}{Print progress information (default: FALSE)}
  \item{parallel}{Use parallel processing (default: TRUE)}
}
\value{
  An S3 object of class \code{"sce"} containing the ensemble model.
}
\examples{
\donttest{
  # Load example data
  data(streamflow_training_10var)
  data(streamflow_testing_10var)

  # Define variables
  Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
  Predictants <- c("Flow")

  # Build SCE model
  sce_model <- sce(
    training_data = streamflow_training_10var,
    x = Predictors,
    y = Predictants,
    mfeature = round(0.5 * length(Predictors)),
    nmin = 5,
    ntree = 48,
    alpha = 0.05,
    resolution = 1000,
    parallel = FALSE
  )

  # Use S3 methods
  print(sce_model)
  summary(sce_model)
  sce_predictions <- predict(sce_model, streamflow_testing_10var)
  sce_importance <- importance(sce_model)
  sce_evaluation <- evaluate(sce_model, streamflow_testing_10var, streamflow_training_10var)
}
}
\seealso{
  \code{\link{sca}}, \code{\link{predict}}, \code{\link{importance}}, \code{\link{evaluate}}
}
