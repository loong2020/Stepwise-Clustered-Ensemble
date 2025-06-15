\name{Streamflow_testing_22var}
\alias{Streamflow_testing_22var}

\title{
  Streamflow Testing Dataset with 22 Variables
}

\description{
  A dataset containing streamflow and related environmental variables for testing purposes.
  This dataset includes additional climate indices and their lagged values.
}

\usage{
data("Streamflow_testing_22var")
}

\format{
  A data frame with the following variables:
  \describe{
    \item{Flow}{Streamflow measurements}
    \item{IPO}{Interdecadal Pacific Oscillation}
    \item{IPO_lag1}{IPO with 1-month lag}
    \item{IPO_lag2}{IPO with 2-month lag}
    \item{Nino3.4}{Nino 3.4 index}
    \item{Nino3.4_lag1}{Nino 3.4 with 1-month lag}
    \item{Nino3.4_lag2}{Nino 3.4 with 2-month lag}
    \item{PDO}{Pacific Decadal Oscillation}
    \item{PDO_lag1}{PDO with 1-month lag}
    \item{PDO_lag2}{PDO with 2-month lag}
    \item{PNA}{Pacific North American pattern}
    \item{PNA_lag1}{PNA with 1-month lag}
    \item{PNA_lag2}{PNA with 2-month lag}
    \item{Precipitation}{Monthly precipitation}
    \item{Precipitation_2Mon}{2-month precipitation}
    \item{Radiation}{Solar radiation}
    \item{Radiation_2Mon}{2-month solar radiation}
    \item{Tmax}{Maximum temperature}
    \item{Tmax_2Mon}{2-month maximum temperature}
    \item{Tmin}{Minimum temperature}
    \item{Tmin_2Mon}{2-month minimum temperature}
    \item{VP}{Vapor pressure}
    \item{VP_2Mon}{2-month vapor pressure}
  }
}

\source{
  The data was collected from environmental monitoring stations and climate indices databases.
}

\examples{
  data("Streamflow_testing_22var")
  str(Streamflow_testing_22var)
} 