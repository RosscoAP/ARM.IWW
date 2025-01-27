# What directory contains the observed CDC data?
CDCdir <- "C:/Users/Pidoto/Desktop/ARM-R-Library/CDC-Input/"
# Set the working directory. This is where the data will be saved to
setwd("C:/Users/Pidoto/Desktop/ARM-R-Library/Output/")

# The following two lines should be left commented out if you are testing locally using Load_Src_Locally.R
library(data.table) # required to be loaded for ARM.IWW
library(ARM.IWW) # load the package - not needed if running functions locally for development

# If also wanting to insert small events, the following three lines are required
library(Rcpp)
library(devtools)
sourceCpp("C:/Users/pidoto/Desktop/ARM-R-Library/ARM.IWW/Small_Events_Rcpp.cpp")

# Collection of stations around Hannover (station 2014)
Stations <- as.character(c(294,2011,2014,5382))

# What parameters are currently set for the model to use?
get_ARM_pars()
# This is how we can set the parameters of the model
set_ARM_pars(Aggregation = 60, WSAmin = 1, DSDmin = 240)
# What model type are we using? Returns either "Callau" or "Pidoto"
get_model_type()

# Loop around stations and perform event definition
for (Station in Stations) {
  # Load observations from CDC source
  Obs <- load_CDC(Station, sourceDir = CDCdir, includeRecent = F)
  # Perform the event definition from observations
  Events <- get_events(Obs, Station)
  # Also determine small events
  Events <- get_events(Obs, Station, SmallEventsMode = T)
}

# First run the model in the Callau mode
set_model_type("Callau")
for (Station in Stations) {
  Marginals <- fit_marginals(Station) # Returns a list of distribution parameters
  Copulas <- fit_copulas(Station) # Returns a list of copula parameters
}
# As the Callau uses a regional empirical copula for WSA-WSD,
# the individual copulas of all stations need to be defined before synthesis
for (Station in Stations) {
  External <- gen_external(Station, Years = 100, RegEmpCopula = Stations) # returns the external time series
  Internal <- gen_internal(Station, Years = 100) # returns the internal time series in units of 1/100th mm
  SmallEvents <- insert_small_events(Station, Years = 100) # returns the internal time series with small events added in units of 1/100th mm
}

# Now run the model in the Pidoto mode
set_model_type("Pidoto")
for (Station in Stations) {
  Marginals <- fit_marginals(Station) # Returns a list of distribution parameters
  Copulas <- fit_copulas(Station) # Returns a list of copula parameters
  External <- gen_external(Station, Years = 100) # returns the external time series
  Internal <- gen_internal(Station, Years = 100) # returns the internal time series in units of 1/100th mm
  SmallEvents <- insert_small_events(Station, Years = 100) # returns the internal time series with small events added in units of 1/100th mm
}

# It is also possible to create realisations of simulations, by including the realisation argument
External <- gen_external("2014", Years = 100, Realisation = 1) # returns the external time series
Internal <- gen_internal("2014", Years = 100, Realisation = 1) # returns the internal time series in units of 1/100th mm
SmallEvents <- insert_small_events("2014", Years = 100, Realisation = 1) # returns the internal time series with small events added in units of 1/100th mm

# It is also possible to override the current model type (see get_model_type()) by using the ModelType argument
# Applies to the functions fit_copulas(), fit_marginals(), gen_external() and gen_internal()
Marginals <- fit_marginals("2014", ModelType = "Callau") # Example
