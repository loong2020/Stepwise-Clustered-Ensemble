\name{Streamflow_datasets}
\alias{Streamflow_training_10var}
\alias{Streamflow_training_22var}
\alias{Streamflow_testing_10var}
\alias{Streamflow_testing_22var}

\title{
  Streamflow Dataset
}

\description{
  These datasets contain streamflow and related environmental variables for training and testing purposes. They are used in examples to demonstrate the SCE package functionality with different levels of complexity.
}

\usage{
data("Streamflow_training_10var")
data("Streamflow_training_22var")
data("Streamflow_testing_10var")
data("Streamflow_testing_22var")
}

\format{
  \strong{Streamflow_training_10var:} Basic environmental variables (12 columns):
  \describe{
    \item{Date}{Date and time of measurement}
    \item{Prcp}{Monthly mean daily precipitation (mm)}
    \item{SRad}{Monthly mean daily solar radiation (W/m^2)}
    \item{Tmax}{Monthly mean daily maximum temperature (°C)}
    \item{Tmin}{Monthly mean daily minimum temperature (°C)}
    \item{VP}{Monthly mean daily vapor pressure (Pa)}
    \item{smlt}{Monthly snowmelt (m)}
    \item{swvl1}{Soil water content layer 1 (m^3/m^3)}
    \item{swvl2}{Soil water content layer 2 (m^3/m^3)}
    \item{swvl3}{Soil water content layer 3 (m^3/m^3)}
    \item{swvl4}{Soil water content layer 4 (m^3/m^3)}
    \item{Flow}{Monthly mean daily streamflow (cfs)}
  }
  
  \strong{Streamflow_training_22var:} Extended variables with climate indices (24 columns):
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
  
  \strong{Testing datasets:} Same structure as corresponding training datasets.
}

\details{
  \strong{Dataset Structure:}
  \itemize{
    \item \strong{10var datasets}: Basic environmental variables (12 columns)
    \item \strong{22var datasets}: Extended variables with climate indices (24 columns)
    \item \strong{Training datasets}: Used for model building
    \item \strong{Testing datasets}: Used for model evaluation
  }
  
  \strong{Climate Indices:} IPO (Interdecadal Pacific Oscillation), Nino3.4 (El Niño), PDO (Pacific Decadal Oscillation), PNA (Pacific North American pattern)
  
  \strong{Data Sources:} ERA5 Land, Daymet, USGS, and climate indices databases
}

\source{
  Environmental monitoring stations, climate indices databases, ERA5 Land, Daymet, and USGS
}
