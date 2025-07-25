\name{RFE_SCE}
\alias{RFE_SCE}
\title{
  Recursive Feature Elimination for SCE Models
}
\description{
  Recursive Feature Elimination for SCE models to identify the most important predictors.
}
\usage{
RFE_SCE(Training_data, Testing_data, Predictors, Predictant, Nmin, Ntree, 
        alpha = 0.05, resolution = 1000, step = 1, verbose = TRUE, 
        parallel = TRUE)
}
\arguments{
  \item{Training_data}{Training dataset}
  \item{Testing_data}{Testing dataset}
  \item{Predictors}{Character vector of predictor names}
  \item{Predictant}{Character vector of predictant names}
  \item{Nmin}{Minimum samples per node}
  \item{Ntree}{Number of trees}
  \item{alpha}{Significance level (default: 0.05)}
  \item{resolution}{Resolution for splitting (default: 1000)}
  \item{step}{Number of predictors to remove per iteration (default: 1)}
  \item{verbose}{Print progress (default: TRUE)}
  \item{parallel}{Use parallel processing (default: TRUE)}
}
\value{
  RFE results with performance metrics and importance scores.
}

% \examples{
% \donttest{
%   # Load example data
%   data(Streamflow_training_22var)
%   data(Streamflow_testing_22var)
%   
%   # Define variables
%   Predictors <- c("IPO","IPO_lag1","IPO_lag2","Nino3.4","Nino3.4_lag1","Nino3.4_lag2",
%                   "PDO","PDO_lag1","PDO_lag2","PNA","PNA_lag1","PNA_lag2",
%                   "Precipitation","Precipitation_2Mon","Radiation","Radiation_2Mon",
%                   "Tmax","Tmax_2Mon","Tmin","Tmin_2Mon","VP","VP_2Mon")
%   Predictants <- c("Flow")
%   
%   # Perform RFE
%   rfe_result <- RFE_SCE(
%     Training_data = Streamflow_training_22var,
%     Testing_data = Streamflow_testing_22var,
%     Predictors = Predictors,
%     Predictant = Predictants,
%     Nmin = 5,
%     Ntree = 48,
%     alpha = 0.05,
%     resolution = 1000,
%     step = 1,
%     verbose = TRUE,
%     parallel = FALSE
%   )
%   
%   # Plot RFE results
%   Plot_RFE(rfe_result)
% }
% }

\seealso{
  \code{\link{Plot_RFE}}, \code{\link{SCE}}, \code{\link{importance}}
} 