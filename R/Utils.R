# to do: set attributes to R objects with DSDmin, WSAmin, DateRange etc.

# defines variables used across functions
arm <- new.env(parent = emptyenv())
arm$DSDmin <- 240
arm$Aggregation <- 60
arm$WSAmin <- 1
arm$ModelType <- "Pidoto"

# Convenience functions for model type
isPidoto <- function() if (arm$ModelType=="Pidoto") return(TRUE) else return(FALSE)
isCallau <- function() if (arm$ModelType=="Callau") return(TRUE) else return(FALSE)

# Convenience text formatting functions
armParText <- function() return(paste0(arm$Aggregation,"min.","DSD",arm$DSDmin,".WSA",arm$WSAmin))
armParTextExt <- function() return(paste0(get_model_type(),".",arm$Aggregation,"min.","DSD",arm$DSDmin,".WSA",arm$WSAmin)) # extended version, includes model type
seasonText <- function(Season) if (Season=="S") return("Summer") else return("Winter")

# Convenience function to pad zeroes
pad0s <- function(X, width=2) return(formatC(X,width=width,format="d",flag="0"))
# Function to load an R object from disk to a given object explicitly
loadRObject <- function(fileLocation) return(eval(parse(text=load(file=fileLocation))))
# Returns either "S" or "W" depending on month
returnSeason <- function(d) return(ifelse(lubridate::month(d) %in% 4:9,"S","W"))
