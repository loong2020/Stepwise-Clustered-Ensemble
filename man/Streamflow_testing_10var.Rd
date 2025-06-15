\name{Streamflow_testing_10var}
\alias{Streamflow_testing_10var}

\title{
  Streamflow Testing Dataset with 10 Variables
}

\description{
  A dataset containing streamflow and related environmental variables for testing purposes.
  This dataset is used in the examples to demonstrate the SCE package functionality.
}

\usage{
data("Streamflow_testing_10var")
}

\format{
  A data frame with the following variables:
  \describe{
    \item{Flow}{Streamflow measurements}
    \item{Prcp}{Precipitation}
    \item{SRad}{Solar radiation}
    \item{Tmax}{Maximum temperature}
    \item{Tmin}{Minimum temperature}
    \item{VP}{Vapor pressure}
    \item{X}{Index variable}
    \item{smlt}{Snow melt}
    \item{swvl1}{Soil water volume layer 1}
    \item{swvl2}{Soil water volume layer 2}
    \item{swvl3}{Soil water volume layer 3}
    \item{swvl4}{Soil water volume layer 4}
  }
}

\source{
  The data was collected from environmental monitoring stations.
}

\examples{
  data("Streamflow_testing_10var")
  str(Streamflow_testing_10var)
}
