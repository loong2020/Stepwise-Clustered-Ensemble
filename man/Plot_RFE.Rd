\name{Plot_RFE}
\alias{Plot_RFE}
\title{
  Plot Recursive Feature Elimination Results
}
\description{
  Plot Recursive Feature Elimination results.
}
\usage{
Plot_RFE(rfe_result, 
         main = "OOB Validation and Testing R2 vs Number of Predictors", 
         col_validation = "blue", 
         col_testing = "red", 
         pch = 16, 
         lwd = 2, 
         cex = 1.2, 
         legend_pos = "bottomleft", 
         ...)
}
\arguments{
  \item{rfe_result}{Result object from RFE_SCE function}
  \item{main}{Plot title}
  \item{col_validation}{Color for validation line}
  \item{col_testing}{Color for testing line}
  \item{pch}{Point character}
  \item{lwd}{Line width}
  \item{cex}{Point size}
  \item{legend_pos}{Legend position}
  \item{...}{Additional arguments}
}
\value{
  Plot showing validation and testing R2 vs number of predictors.
}

\seealso{
  \code{\link{RFE_SCE}}
} 