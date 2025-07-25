\name{Air_quality_datasets}
\alias{Air_quality_training}
\alias{Air_quality_testing}

\title{
  Air Quality Dataset
}

\description{
  These datasets contain air quality measurements for training and testing purposes. They include various air pollutant concentrations and meteorological variables measured at different locations and times.
}

\usage{
data("Air_quality_training")
data("Air_quality_testing")
}

\format{
  Both datasets are data frames with 8760 rows and 12 variables:
  \describe{
    \item{Date}{Date and time of measurement (POSIXct format)}
    \item{PM2.5}{Particulate matter with diameter less than 2.5 micrometers (\\mu g/m^3)}
    \item{PM10}{Particulate matter with diameter less than 10 micrometers (\\mu g/m^3)}
    \item{SO2}{Sulfur dioxide concentration (\\mu g/m^3)}
    \item{NO2}{Nitrogen dioxide concentration (\\mu g/m^3)}
    \item{CO}{Carbon monoxide concentration (\\mu g/m^3)}
    \item{O3}{Ozone concentration (\\mu g/m^3)}
    \item{TEMP}{Temperature (\\textdegree C)}
    \item{PRES}{Atmospheric pressure (hPa)}
    \item{DEWP}{Dew point temperature (\\textdegree C)}
    \item{RAIN}{Precipitation amount (mm)}
    \item{WSPM}{Wind speed (m/s)}
  }
}

\details{
  \strong{Dataset Differences:}
  \itemize{
    \item \code{Air_quality_training}: Used for training SCA and SCE models
    \item \code{Air_quality_testing}: Used for testing trained models
  }
  
  \strong{Variable Descriptions:}
  \itemize{
    \item \strong{PM2.5, PM10}: Particulate matter concentrations, important indicators of air quality
    \item \strong{SO2, NO2, CO, O3}: Major air pollutants regulated by environmental agencies
    \item \strong{TEMP, PRES, DEWP}: Meteorological variables affecting air quality
    \item \strong{RAIN, WSPM}: Weather conditions that influence pollutant dispersion
  }
}

\source{
  Air quality monitoring stations
}