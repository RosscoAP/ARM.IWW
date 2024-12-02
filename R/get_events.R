###############################################################################
#                        FUNCTIONS USED IN THIS SCRIPT                        #
###############################################################################

# Function for determining the external structure of the time series based on wet and dry periods
examineExternalStructure <- function(sourceTS) {
  extStructure <- rle(sourceTS$Rain) # run length encoding function - compute the lengths and values of runs of equal values in a vector
  extStructure <- data.frame(Length=extStructure[[1]],Rain=extStructure[[2]]) # convert to data frame
  extStructure$endIndex <- cumsum(extStructure$Length)
  extStructure$startIndex <- extStructure$endIndex - extStructure$Length + 1
  extStructure$NAs <- unlist(mapply(eventNAs, extStructure$startIndex, extStructure$endIndex, MoreArgs=list(sourceTS=sourceTS)))
  extStructure$WSA[extStructure$Rain] <- unlist(mapply(eventWSA, extStructure$startIndex, extStructure$endIndex, MoreArgs=list(sourceTS=sourceTS)))[extStructure$Rain]
  extStructure$WSP[extStructure$Rain] <- unlist(mapply(eventWSP, extStructure$startIndex, extStructure$endIndex, MoreArgs=list(sourceTS=sourceTS)))[extStructure$Rain]
  extStructure$WSPT[extStructure$Rain] <- unlist(mapply(eventTimeToPeak, extStructure$startIndex, extStructure$endIndex, extStructure$WSP, MoreArgs=list(sourceTS=sourceTS)))[extStructure$Rain]*arm$Aggregation-arm$Aggregation # in minutes
  eventsPcpWSP <<- unlist(mapply(eventPcpByTimeStepWSP, extStructure$startIndex[extStructure$Rain], extStructure$endIndex[extStructure$Rain], MoreArgs=list(sourceTS=sourceTS)))
  eventsPcpWSA <<- unlist(mapply(eventPcpByTimeStepWSA, extStructure$startIndex[extStructure$Rain], extStructure$endIndex[extStructure$Rain], MoreArgs=list(sourceTS=sourceTS)))
  return(extStructure)
}
# Functions to calculate the WSA, WSP, TimeToPeak for an event
eventNAs <- function(startIndex,endIndex,sourceTS) return(sum(sourceTS$isNA[startIndex:endIndex]))
eventWSA <- function(startIndex,endIndex,sourceTS) return(sum(sourceTS$Pcp[startIndex:endIndex], na.rm = TRUE))
eventWSP <- function(startIndex,endIndex,sourceTS) {
  if (sum(sourceTS$Pcp[startIndex:endIndex], na.rm = TRUE)==0) return(NA) else return(max(sourceTS$Pcp[startIndex:endIndex], na.rm = TRUE))
}
eventPcpByTimeStepWSP <- function(startIndex,endIndex,sourceTS) return( sourceTS$Pcp[startIndex:endIndex] / max(sourceTS$Pcp[startIndex:endIndex], na.rm = TRUE) )
eventPcpByTimeStepWSA <- function(startIndex,endIndex,sourceTS) return( sourceTS$Pcp[startIndex:endIndex] / sum(sourceTS$Pcp[startIndex:endIndex], na.rm = TRUE) )
eventTimeToPeak <- function(startIndex,endIndex,WSP,sourceTS) {
  #print(paste(startIndex,endIndex,WSP))
  if (is.na(WSP)) return(NA) else return(sample(rep(which(!is.na(match(sourceTS$Pcp[startIndex:endIndex],WSP))),2),1))
  # the rep() part ensures that the vector to sample from is never of length 1, as this can lead to undesirable results (see documentation for sample() )
}

###############################################################################
#                             EVENT SELECTION CODE                            #
###############################################################################

#' Title
#'
#' @param PcpTS
#' @param StationID
#' @param SmallEventsMode
#' @param WSImin
#' @param DSDMax
#' @param saveRObjects
#' @param saveText
#' @param includeEventTS
#'
#' @return
#' @export
#'
#' @examples
get_events <- function(PcpTS, StationID,
                       SmallEventsMode = FALSE,
                       WSImin = 0.0099, # [mm/timestep]
                       DSDMax = 6*7*24*60, # [min] Used to kick out DSD outliers (replaced by NAs, previously the median DSD)
                       saveRObjects = TRUE,
                       saveText = FALSE,
                       includeEventTS = FALSE
                       ) {

  # How many timesteps account for arm$DSDmin?
  DSD_MinSteps <- arm$DSDmin/arm$Aggregation #(min/min)
  # Check that DSD_MinSteps is an integer >= 1
  if ((DSD_MinSteps<1)|(as.integer(DSD_MinSteps)!=DSD_MinSteps)) stop("Please check function arguments. arm$DSDmin is not a multiple of timestep.")

  ############################
  # Data preparation / cleanup
  totalRainfall <- sum(PcpTS$Pcp, na.rm = TRUE)
  print.noquote(paste("## Event definition:", StationID, "##"))
  print.noquote(paste("Total Rainfall:", round(totalRainfall,2),"mm"))
  print.noquote(paste("Rainfall <= WSImin:", round(sum(PcpTS$Pcp[PcpTS$Pcp<=WSImin], na.rm = TRUE),2),"mm"))
  # If there is no rainfall in the time series, abort the script
  try(if(totalRainfall==0) stop("Time series includes no rainfall events."))
  # Correcting the time series from the time step with intensity>WSImin.
  PcpTS$Pcp[round(PcpTS$Pcp,2) <= WSImin] <- 0

  ##################
  # DSDmin condition
  # return the external structure
  events <- examineExternalStructure(PcpTS)
  # First find all the dry periods that don't meet the arm$DSDmin condition
  shortDSD <- events[(!events$Rain & events$Length<DSD_MinSteps),]
  # Set these too short dry periods periods to wet periods
  PcpTS$Rain[unlist(mapply(":", shortDSD$startIndex, shortDSD$endIndex))] <- TRUE
  # Re-run the external structre
  events <- examineExternalStructure(PcpTS)

  ##################
  # arm$WSAmin condition
  if (SmallEventsMode) {
    # First find all the wet periods that are above the arm$WSAmin condition
    highWSA <- events[(events$Rain & events$WSA>=arm$WSAmin),]
    # Set these too high wet periods periods to dry periods
    PcpTS$Rain[unlist(mapply(":", highWSA$startIndex, highWSA$endIndex))] <- FALSE
    # For small events, we ignore the effect of arm$DSDmin
    PcpTS$Rain[PcpTS$Pcp==0] <- FALSE # set timesteps with zero rainfall to dry
    # Reverse the periods that were set to wet above due to being under arm$DSDmin
    PcpTS$Rain[unlist(mapply(":", shortDSD$startIndex, shortDSD$endIndex))] <- FALSE
  } else {
    # First find all the wet periods that don't meet the arm$WSAmin condition
    lowWSA <- events[(events$Rain & events$WSA<arm$WSAmin),]
    # Set these too low wet periods periods to dry periods
    PcpTS$Rain[unlist(mapply(":", lowWSA$startIndex, lowWSA$endIndex))] <- FALSE
  }

  # Re-run the external structre
  events <- examineExternalStructure(PcpTS)
  print.noquote(paste("Rainfall >= WSAmin:", round(sum(events$WSA, na.rm = TRUE),2),"mm"))
  print.noquote(paste("% of total rainfall:", round(sum(events$WSA, na.rm = TRUE)/totalRainfall,2)*100,"%"))

  ####################
  # Dataframe clean-up
  # Remove the very first event, as we are not sure as to it's true length
  events <- events[2:nrow(events),]
  # Then for consistency, we want the table to be in the order WET-DRY
  if(!events$Rain[1]) events <- events[2:nrow(events),]
  # At the same time, make sure the last entry is a dry event
  if(events$Rain[nrow(events)]) events <- events[1:nrow(events)-1,]
  # Move the DSD information up to the wet spell part
  events$DSD <- NA
  events$DSD[seq(1,nrow(events),by=2)] <- events$Length[seq(2,nrow(events),by=2)]*arm$Aggregation # in minutes
  # Same for dry NA count
  events$DryNAs <- NA; events$DryNAs[seq(1,nrow(events),by=2)] <- events$NAs[seq(2,nrow(events),by=2)]
  # Shift wet NA count to end of data.frame
  events$WetNAs <- events$NAs; events$NAs <- NULL
  # Now remove every second row (the dry part)
  events <- events[seq(1,nrow(events),by=2),]
  events$Rain <- NULL
  events$WSD <- events$Length*arm$Aggregation # in minutes
  events$Length <- NULL
  events$StartDateTime <- PcpTS$DateTime[events$startIndex] # Obtain the starting dates
  events$WSI <- round(events$WSA / events$WSD,digits=5)
  events$TS_Index <- events$startIndex
  if (includeEventTS) events$Pcp <- lapply(X = mapply(":", events$startIndex, events$endIndex), FUN = function(x) { PcpTS$Pcp[x] })
  events$endIndex <- NULL
  events$startIndex <- NULL
  events$Season <- "S" # Create a season variable - either "S" or "W"
  events$Season[which(is.na(match(lubridate::month(events$StartDateTime), 4:9 )))] <- "W" # assign to winter if between April and September
  events$Season <- as.factor(events$Season)
  # replace DSD outliers
  print.noquote(paste("DSD outliers:", sum(events$DSD>DSDMax)))
  events$DSD[which(events$DSD>DSDMax)] <- NA # median(events$DSD)
  print.noquote(paste("Rainfall events:", nrow(events)))
  #print.noquote(lubridate::year(events$StartDateTime[1]) + lubridate::month(events$StartDateTime[1])/12) # number of years

  ################
  # Output results
  # create a sub-folder within the output folder made up of definition parameters
  outputSubfolder <- paste0(armParText(),ifelse(SmallEventsMode,"_SmallEvents",""),"/")
  if (saveRObjects | saveText) if (!dir.exists(paste0("Events/",outputSubfolder))) {
    dir.create("Events/")
    dir.create(paste0("Events/",outputSubfolder))
  }
  # Save CSV text object if chosen
  if (saveText) {
    if ( includeEventTS) write.csv(events[,-c("Pcp")], file=paste0("Events/",outputSubfolder,StationID,".txt"), row.names = FALSE)
    if (!includeEventTS) write.csv(events[], file=paste0("Events/",outputSubfolder,StationID,".txt"), row.names = FALSE)
  }
  # save RObjects if chosen
  if (saveRObjects) save(events, file=paste0("Events/",outputSubfolder,StationID,".RData"))

  return(events)

}
