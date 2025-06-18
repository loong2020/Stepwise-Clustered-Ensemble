\name{Air_quality_training}
\alias{Air_quality_training}

\title{
	Air Quality Training Dataset
}

\description{
	This dataset contains air quality measurements for training purposes. It includes various air pollutant concentrations measured at different locations and times.
}

\usage{
data("Air_quality_training")
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
data(Air_quality_training)
head(Air_quality_training)
}