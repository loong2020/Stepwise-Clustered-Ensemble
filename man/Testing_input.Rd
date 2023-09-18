\name{Testing_input}
\alias{Testing_input}
\alias{Testing_input}

\title{
	Testing dataset for Example
}

\description{
	This testing dataset is used as an example to illustrate the SCE model.
}

\usage{
data("Testing_input")
}
\format{
  A data frame with X rows and Y columns:
  \describe{
    \item{Date}{Time.}
    \item{Prcp}{Monthly mean daily precipitation (mm), derived from DAYMET.}
    \item{SRad}{Monthly mean daily short-wave solar radiation (W/m^2), derived from DAYMET.}
    \item{Tmax}{Monthly mean daily maximal temperature (°C), derived from DAYMET.}
    \item{Tmin}{Monthly mean daily minimal temperature (°C), derived from DAYMET.}
    \item{VP}{Monthly mean daily vapor pressure (Pa), derived from DAYMET.}
    \item{smlt}{Monthly snowmelt sum (m), derived from ERA5.}
    \item{swvl1}{Monthly volumetric soil water layer 1 (m³/m³), derived from ERA5.}
    \item{swvl2}{Monthly volumetric soil water layer 2 (m³/m³), derived from ERA5.}
    \item{swvl3}{Monthly volumetric soil water layer 3 (m³/m³), derived from ERA5.}
    \item{swvl4}{Monthly volumetric soil water layer 4 (m³/m³), derived from ERA5.}
    \item{Flow}{Monthly mean daily streamflow rate (cubic feet per second), derived from USGS.}
  }
}
\source{
  The data is sourced from ERA5, DAYMET, and USGS.
}
\examples{
  data(Testing_input)
  head(Testing_input)
}
