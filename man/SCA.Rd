\name{SCA}
\alias{SCA}
\title{
  Stepwise Cluster Analysis (SCA)
}
\description{
  Builds a single Stepwise Cluster Analysis (SCA) tree model that recursively partitions the data space based on Wilks' Lambda statistic.
}
\usage{
SCA(Training_data, X, Y, Nmin, alpha = 0.05, resolution = 1000, verbose = FALSE)
}
\arguments{
  \item{Training_data}{A data.frame containing the training data}
  \item{X}{Character vector of predictor variable names}
  \item{Y}{Character vector of predictant variable names}
  \item{Nmin}{Minimum number of samples in a leaf node}
  \item{alpha}{Significance level for clustering (default: 0.05)}
  \item{resolution}{Resolution for splitting (default: 1000)}
  \item{verbose}{Print progress information (default: FALSE)}
}
\value{
  An S3 object of class "SCA" containing the tree model.
}
\examples{
\donttest{
  # Load example data
  data(Streamflow_training_10var)
  data(Streamflow_testing_10var)
  
  # Define variables
  Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
  Predictants <- c("Flow")
  
  # Build SCA model
  sca_model <- SCA(
    Training_data = Streamflow_training_10var,
    X = Predictors,
    Y = Predictants,
    Nmin = 5,
    alpha = 0.05,
    resolution = 1000
  )
  
  # Use S3 methods
  print(sca_model)
  summary(sca_model)
  sca_predictions <- predict(sca_model, Streamflow_testing_10var)
  sca_importance <- importance(sca_model)
  sca_evaluation <- evaluate(sca_model, Streamflow_testing_10var, Streamflow_training_10var)
}
}
\seealso{
  \code{\link{SCE}}, \code{\link{predict}}, \code{\link{importance}}, \code{\link{evaluate}}
} 