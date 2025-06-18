\name{Air_quality_testing}
\alias{Air_quality_testing}

\title{
  Example Air Quality Testing Dataset for the SCE Model
}

\description{
  The "Air_quality_testing" dataset serves as an example dataset to illustrate the functioning and application of the SCE model within the SCE package. It contains various air quality parameters and meteorological variables measured on a monthly scale, which are used as inputs in the SCE model for model testing.
}

\usage{
data("Air_quality_testing")
}

\format{
	A data frame with 8760 rows and 7 variables:
	\describe{
		\item{Date}{Date and time of measurement (POSIXct format)}
		\item{PM2.5}{Particulate matter with diameter less than 2.5 micrometers (\\mu g/m^3).}
		\item{PM10}{Particulate matter with diameter less than 10 micrometers (\\mu g/m^3).}
		\item{SO2}{Sulfur dioxide concentration (\\mu g/m^3).}
		\item{NO2}{Nitrogen dioxide concentration (\\mu g/m^3).}
		\item{CO}{Carbon monoxide concentration (\\mu g/m^3).}
		\item{O3}{Ozone concentration (\\mu g/m^3).}
		\item{TEMP}{Temperature (\\textdegree C).}
		\item{PRES}{Atmospheric pressure (hPa).}
		\item{DEWP}{Dew point temperature (\\textdegree C).}
		\item{RAIN}{Precipitation amount (mm).}
		\item{WSPM}{Wind speed (m/s).}
	}
}

\source{
	Air quality monitoring stations
}

\examples{
data(Air_quality_testing)
head(Air_quality_testing)
}