## Used for local testing without needing to build and load library

# Where are the source files located?
srcDir <- "C:/Users/pidoto/Desktop/ARM-R-Library/ARM.IWW/R/" # must include trailing forward slash

# Now all the relevant scripts in that directory will be loaded into the global environment
source(paste0(srcDir,"Utils.R"))
source(paste0(srcDir,"gen_internal.R"))
source(paste0(srcDir,"get_events.R"))
source(paste0(srcDir,"load_CDC.R"))
source(paste0(srcDir,"gen_external.R"))
source(paste0(srcDir,"fit_copulas.R"))
source(paste0(srcDir,"fit_marginals.R"))
source(paste0(srcDir,"model_type.R"))
source(paste0(srcDir,"ARM_params.R"))
source(paste0(srcDir,"season.R"))

# Load required libraries
library(data.table)
library(lmomco)
library(copula)
library(lubridate)
