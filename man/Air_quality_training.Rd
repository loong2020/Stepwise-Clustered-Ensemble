\name{Air_quality_training}
\alias{Air_quality_training}
\alias{Air_quality_training}

\title{
  Example Air Quality Training Dataset for the SCE Model
}

\description{
  The "Air_quality_training" dataset serves as an example dataset to illustrate the functioning and application of the SCE model within the SCE package. It contains various air quality parameters and meteorological variables measured on a monthly scale, which are used as inputs in the SCE model for model training.
}

\usage{
data("Air_quality_training")
}

\format{
  A matrix containing the following variables:

  \describe{
    \item{PM2.5}{Particulate matter with diameter less than 2.5 micrometers (μg/m³).}
    \item{PM10}{Particulate matter with diameter less than 10 micrometers (μg/m³).}
    \item{SO2}{Sulfur dioxide concentration (μg/m³).}
    \item{NO2}{Nitrogen dioxide concentration (μg/m³).}
    \item{CO}{Carbon monoxide concentration (μg/m³).}
    \item{O3}{Ozone concentration (μg/m³).}
    \item{TEMP}{Temperature (°C).}
    \item{PRES}{Atmospheric pressure (hPa).}
    \item{DEWP}{Dew point temperature (°C).}
    \item{RAIN}{Precipitation amount (mm).}
    \item{WSPM}{Wind speed (m/s).}
  }
}

\source{
  The data is compiled from air quality monitoring station and meteorological station at Changping District, Beijing, providing comprehensive measurements of various air pollutants and weather conditions.
}

\examples{
  # Load the "Air_quality_training" dataset and display the first few rows
  data(Air_quality_training)
  head(Air_quality_training)
} 