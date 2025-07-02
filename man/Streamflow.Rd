\name{Streamflow_training_10var}
\alias{Streamflow_training_10var}
\alias{Streamflow_training_22var}
\alias{Streamflow_testing_10var}
\alias{Streamflow_testing_22var}

\title{
  Streamflow Datasets
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
  \strong{Streamflow_training_10var:} A data frame with basic environmental variables:
  \describe{
    \item{Date}{The date and time of the data point}
    \item{Prcp}{The monthly mean daily precipitation measured in millimeters (mm), derived from the Daymet dataset}
    \item{SRad}{The monthly mean daily short-wave solar radiation measured in Watts per square meter (W/m^2), sourced from the Daymet dataset}
    \item{Tmax}{The monthly mean daily maximal temperature recorded in degrees Celsius, taken from the Daymet dataset}
    \item{Tmin}{The monthly mean daily minimal temperature recorded in degrees Celsius, also derived from the Daymet dataset}
    \item{VP}{The monthly mean daily vapor pressure measured in Pascals (Pa), obtained from the Daymet dataset}
    \item{smlt}{The sum of monthly snowmelt measurements in meters (m), taken from the ERA5 land dataset}
    \item{swvl1}{The volumetric soil water content in layer 1 measured in cubic meters per cubic meter (m^3/m^3), sourced from the ERA5 land dataset}
    \item{swvl2}{The volumetric soil water content in layer 2, measured similarly to swvl1, sourced from the ERA5 land dataset}
    \item{swvl3}{The volumetric soil water content in layer 3, measured similarly to swvl1, sourced from the ERA5 land dataset}
    \item{swvl4}{The volumetric soil water content in layer 4, measured similarly to swvl1, sourced from the ERA5 land dataset}
    \item{Flow}{The monthly mean daily streamflow rate measured in cubic feet per second (cfs), provided by the United States Geological Survey (USGS)}
  }
  
  \strong{Streamflow_training_22var:} A data frame with extended variables including climate indices:
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
  
  \strong{Streamflow_testing_10var:} A data frame with basic environmental variables (same structure as training):
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
  
  \strong{Streamflow_testing_22var:} A data frame with extended variables including climate indices (same structure as training):
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

\details{
  \strong{Dataset Categories:}
  \itemize{
    \item \strong{Training Datasets}: Used for building SCA and SCE models
      \itemize{
        \item \code{Streamflow_training_10var}: Basic dataset with 12 variables, suitable for introductory examples
        \item \code{Streamflow_training_22var}: Extended dataset with 24 variables, includes climate indices and lagged values
      }
    \item \strong{Testing Datasets}: Used for evaluating trained models
      \itemize{
        \item \code{Streamflow_testing_10var}: Basic dataset with 12 variables, matches training structure
        \item \code{Streamflow_testing_22var}: Extended dataset with 24 variables, matches training structure
      }
  }
  
  \strong{Variable Categories:}
  \itemize{
    \item \strong{Hydrological}: Flow, Precipitation, Snow melt, Soil water volumes
    \item \strong{Meteorological}: Temperature (max/min), Solar radiation, Vapor pressure
    \item \strong{Climate Indices}: IPO, Nino3.4, PDO, PNA (with lagged versions)
    \item \strong{Time Aggregations}: 2-month averages for key variables
  }
  
  \strong{Climate Indices:}
  \itemize{
    \item \strong{IPO}: Interdecadal Pacific Oscillation - long-term climate pattern
    \item \strong{Nino3.4}: El Niño-Southern Oscillation index
    \item \strong{PDO}: Pacific Decadal Oscillation - long-term ocean temperature pattern
    \item \strong{PNA}: Pacific North American pattern - atmospheric circulation pattern
  }
  
  \strong{Data Sources:}
  The data is compiled from various recognized sources including:
  \itemize{
    \item ERA5 Land: A global land-surface dataset at 9km resolution, available from the Copernicus Climate Change Service
    \item Daymet Version 4: Daily Surface Weather and Climatological Summaries
    \item United States Geological Survey (USGS): A scientific agency of the United States government that studies natural resources, natural hazards, and the landscape of the United States
    \item Climate indices databases for the extended datasets
  }
}

\source{
  Environmental monitoring stations, climate indices databases, ERA5 Land, Daymet, and USGS
}
