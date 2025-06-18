\name{OOB_validation}
\alias{OOB_validation}

\title{
  Calculate Out-of-Bag (OOB) Validation Predictions
}

\description{
  This function calculates the out-of-bag (OOB) validation predictions for an SCE model.
  OOB predictions are made using only the trees that did not use a particular observation
  during training, providing an unbiased estimate of the model's performance.
}

\usage{
OOB_validation(model)
}

\arguments{
  \item{model}{
    A trained SCE model object generated through the \code{SCE} function.
    The model must contain OOB information in its Tree_Info component.
  }
}

\value{
  A data.frame containing the OOB predictions for each predictant.
  The number of rows equals the number of training samples, and the columns
  correspond to the predictant variables.
}

\details{
  The OOB validation process involves:
  \enumerate{
    \item Extracting OOB indices and predictions from each tree
    \item Weighting the predictions based on tree weights
    \item Calculating weighted means for each sample
  }
  
  This function is typically called internally by \code{Model_simulation}
  and is not usually called directly by users.
}

\author{
  Kailong Li <lkl98509509@gmail.com>
}

\seealso{
  \code{\link{SCE}}
} 