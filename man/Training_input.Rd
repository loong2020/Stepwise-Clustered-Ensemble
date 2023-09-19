\name{Training_input}
\alias{Training_input}
\alias{Training_input}

\title{
  Example Training Dataset for the SCE Model
}

\description{
  The "Training_input" dataset serves as an example dataset to illustrate the functioning and application of the SCE model within the SCE package. It contains various environmental variables measured on a monthly scale, which are used as inputs in the SCE model for model training.
}

\usage{
data("Training_input")
}

\format{
  A data frame with a number of rows (X) and columns (Y) containing the following variables:

  \describe{
    \item{Date}{The date and time of the data point.}
    \item{Prcp}{The monthly mean daily precipitation measured in millimeters (mm). This data is derived from the Daymet dataset.}
    \item{SRad}{The monthly mean daily short-wave solar radiation measured in Watts per square meter (W/m^2), sourced from the Daymet dataset.}
    \item{Tmax}{The monthly mean daily maximal temperature recorded in degrees Celsius, taken from the Daymet dataset.}
    \item{Tmin}{The monthly mean daily minimal temperature recorded in degrees Celsius, also derived from the Daymet dataset.}
    \item{VP}{The monthly mean daily vapor pressure measured in Pascals (Pa), obtained from the Daymet dataset.}
    \item{smlt}{The sum of monthly snowmelt measurements in meters (m), taken from the ERA5 land dataset.}
    \item{swvl1}{The volumetric soil water content in layer 1 measured in cubic meters per cubic meter (m^3/m^3), sourced from the ERA5 land dataset.}
    \item{swvl2}{The volumetric soil water content in layer 2, measured similarly to swvl1, sourced from the ERA5 land dataset.}
    \item{swvl3}{The volumetric soil water content in layer 3, measured similarly to swvl1, sourced from the ERA5 land dataset.}
    \item{swvl4}{The volumetric soil water content in layer 4, measured similarly to swvl1, sourced from the ERA5 land dataset.}
    \item{Flow}{The monthly mean daily streamflow rate measured in cubic feet per second (cfs), provided by the United States Geological Survey (USGS).}
  }
}

\source{
  The data is compiled from various recognized sources including:
  \itemize{
    \item ERA5 Land: A global land-surface dataset at 9km resolution, available from the Copernicus Climate Change Service.
    \item Daymet Version 4: Daily Surface Weather and Climatological Summaries.
    \item United States Geological Survey (USGS): A scientific agency of the United States government that studies natural resources, natural hazards, and the landscape of the United States.
  }
}

\examples{
  # Load the "Training_input" dataset and display the first few rows
  data(Training_input)
  head(Training_input)
}
