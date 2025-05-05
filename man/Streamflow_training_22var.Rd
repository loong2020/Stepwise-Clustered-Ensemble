\name{Streamflow_training_22var}
\alias{Streamflow_training_22var}
\alias{Streamflow_training_22var}

\title{
  Example Streamflow Training Dataset for the SCE Model (22 variables)
}

\description{
  The "Streamflow_training_22var" dataset serves as an example dataset to illustrate the functioning and application of the SCE model within the SCE package. It contains monthly mean values of meteorological data from the CAMELS dataset and four commonly used large-scale climatic indices (Nino3.4, PDO, IPO, and PNA). To reflect initial catchment conditions and lagged impacts, the dataset includes both current values and 2-month moving averages of meteorological variables, as well as current and lagged values (up to 2 months) of climatic indices.
}

\usage{
data("Streamflow_training_22var")
}

\format{
  A data frame containing the following variables:

  \describe{
    \item{Date}{The date and time of the data point.}
    \item{Precipitation}{The monthly mean precipitation (mm).}
    \item{Radiation}{The monthly mean solar radiation (W/m^2).}
    \item{Tmax}{The monthly mean maximum temperature (°C).}
    \item{Tmin}{The monthly mean minimum temperature (°C).}
    \item{VP}{The monthly mean vapor pressure (Pa).}
    \item{Precipitation_2Mon}{The 2-month moving average of precipitation (mm).}
    \item{Radiation_2Mon}{The 2-month moving average of solar radiation (W/m^2).}
    \item{Tmax_2Mon}{The 2-month moving average of maximum temperature (°C).}
    \item{Tmin_2Mon}{The 2-month moving average of minimum temperature (°C).}
    \item{VP_2Mon}{The 2-month moving average of vapor pressure (Pa).}
    \item{PNA}{The Pacific North American index at time t.}
    \item{Nino3.4}{The Nino3.4 index at time t.}
    \item{IPO}{The Interdecadal Pacific Oscillation index at time t.}
    \item{PDO}{The Pacific Decadal Oscillation index at time t.}
    \item{PNA_lag1}{The Pacific North American index at time t-1.}
    \item{Nino3.4_lag1}{The Nino3.4 index at time t-1.}
    \item{IPO_lag1}{The Interdecadal Pacific Oscillation index at time t-1.}
    \item{PDO_lag1}{The Pacific Decadal Oscillation index at time t-1.}
    \item{PNA_lag2}{The Pacific North American index at time t-2.}
    \item{Nino3.4_lag2}{The Nino3.4 index at time t-2.}
    \item{IPO_lag2}{The Interdecadal Pacific Oscillation index at time t-2.}
    \item{PDO_lag2}{The Pacific Decadal Oscillation index at time t-2.}
    \item{Flow}{The monthly mean daily streamflow rate measured in cubic feet per second (cfs).}
  }
}

\source{
  The data is compiled from various recognized sources including:
  \itemize{
    \item CAMELS (Catchment Attributes and Meteorology for Large-sample Studies): A large-sample hydrology dataset.
    \item Nino3.4 Index: Based on Trenberth (1997).
    \item Pacific Decadal Oscillation (PDO): Based on Mantua et al. (1997).
    \item Interdecadal Pacific Oscillation (IPO): Based on Mantua et al. (1997).
    \item Pacific North American Index (PNA): Based on Leathers et al. (1991).
  }
}

\references{
  \itemize{
    \item Trenberth, K. E. (1997). The definition of El Niño. Bulletin of the American Meteorological Society, 78(12), 2771-2777.
    \item Mantua, N. J., Hare, S. R., Zhang, Y., Wallace, J. M., & Francis, R. C. (1997). A Pacific interdecadal climate oscillation with impacts on salmon production. Bulletin of the American Meteorological Society, 78(6), 1069-1079.
    \item Leathers, D. J., Yarnal, B., & Palecki, M. A. (1991). The Pacific/North American teleconnection pattern and United States climate. Part I: Regional temperature and precipitation associations. Journal of Climate, 4(5), 517-528.
  }
}

\examples{
  # Load the "Streamflow_training_22var" dataset and display the first few rows
  data(Streamflow_training_22var)
  head(Streamflow_training_22var)
} 