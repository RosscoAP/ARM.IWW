#' Title
#'
#' @param StationID
#' @param ObsDateFrom
#' @param ObsDateTo
#' @param Years
#' @param Realisation
#' @param saveObjects
#' @param ModelType
#'
#' @return
#' @export
#'
#' @examples
insert_small_events <- function(StationID,
                         ObsDateFrom = NULL,
                         ObsDateTo = NULL,
                         Years,
                         Realisation = NA,
                         saveObjects = TRUE,
                         ModelType = get_model_type()
) {

  # To do: check that input files exist

  set_model_type(ModelType)

  print.noquote(paste("## Insert small events:", StationID, "##"))
  startTime <- Sys.time()

  ##########################################
  # First work out small events parameters #
  ##########################################

  # Function to find closest large event from small event and return time to event [minutes]
  FindLargeEvent <- function(Start, End) {
    after <- minPositive(as.integer(Start) - as.integer(LargeEvents$EndDateTime) )/60
    before <- minPositive(as.integer(LargeEvents$StartDateTime) - as.integer(End) )/60
    # closest large event occurs before small event
    if (before<after) return( -before ) else return( after )
  }
  minPositive = function(x) min(x[x > 0]) # return minimum positive value

  # Load events for station
  SmallEvents <- events <- loadRObject(paste0("Events/",armParText(),"_SmallEvents/",StationID,".RData"))
  LargeEvents <- events <- loadRObject(paste0("Events/",armParText(),"/",StationID,".RData"))

  # Crop the time series if desired
  if (!is.null(ObsDateFrom)) { # TODO: check validity of ObsDateFrom
    SmallEvents <- SmallEvents[SmallEvents$StartDateTime >= ObsDateFrom,]
    LargeEvents <- LargeEvents[LargeEvents$StartDateTime >= ObsDateFrom,]
  }
  if (!is.null(ObsDateTo)) { # TODO: check validity of ObsDateTo
    SmallEvents <- SmallEvents[SmallEvents$StartDateTime <= ObsDateTo,]
    LargeEvents <- LargeEvents[LargeEvents$StartDateTime <= ObsDateTo,]
  }

  LargeEvents$Year <- lubridate::year(LargeEvents$StartDateTime)
  SmallEvents$Year <- lubridate::year(SmallEvents$StartDateTime)

  # Calculate end times of both small and large events
  SmallEvents$EndDateTime <- SmallEvents$StartDateTime + SmallEvents$WSD*60
  LargeEvents$EndDateTime <- LargeEvents$StartDateTime + LargeEvents$WSD*60

  # Remove NA values
  SmallEvents <- SmallEvents[!is.na(SmallEvents$EndDateTime),]
  LargeEvents <- LargeEvents[!is.na(LargeEvents$EndDateTime),]

  # Calculate small event distance to next large event
  SmallEvents$Distance <- mapply(FUN = FindLargeEvent, Start = SmallEvents$StartDateTime, End = SmallEvents$EndDateTime)

  # Calculate mean annual large/small event volumes
  ObsYears <- (sum(LargeEvents$DSD, na.rm = T) + sum(LargeEvents$WSD, na.rm = T)) / (24*60*365)

  Small.Summer <- sum(SmallEvents$WSA[SmallEvents$Season=="S"] ) / ObsYears
  Small.Winter <- sum(SmallEvents$WSA[SmallEvents$Season=="W"] ) / ObsYears
  Large.Summer <- sum(LargeEvents$WSA[LargeEvents$Season=="S"] ) / ObsYears
  Large.Winter <- sum(LargeEvents$WSA[LargeEvents$Season=="W"] ) / ObsYears

  # Ratios between large and small events for both summer and winter
  Summer.Ratio <- (Small.Summer / Large.Summer)*100
  Winter.Ratio <- (Small.Winter / Large.Winter)*100
  Summer.Volume <- Small.Summer
  Winter.Volume <- Small.Winter

  # Separate by season
  Summer.SmallEvents <- SmallEvents[SmallEvents$Season=="S",c("WSA","WSD","Distance")]
  Winter.SmallEvents <- SmallEvents[SmallEvents$Season=="W",c("WSA","WSD","Distance")]

  ###################################################################
  # Finished calculation or small event parameters, start inserting #
  ###################################################################

  # move seasons to list for better handling
  SmallEvents <- list()
  SmallEvents[["S"]] <- Summer.SmallEvents; SmallEvents[["W"]] <- Winter.SmallEvents; rm(Winter.SmallEvents,Summer.SmallEvents);

  # Load external and internal time series
  Pcp <- loadRObject(paste0("Internal/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"))$Pcp
  External <- loadRObject(paste0("External/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"))

  # cut down to important vars
  External <- External[,c("Date","WSD","WSA","DSD")]
  External$Season <- returnSeason(External$Date)
  # convert to timesteps
  External[,c("WSD","DSD")] <- External[,c("WSD","DSD")]/arm$Aggregation
  # calculate the starting timestep of dry steps
  External$StartDSD <- cumsum(External$WSD) + cumsum(External$DSD) - External$DSD # no need to add 1 as indexing in C++ starts at 0
  # Indexes of events for summer and winter
  EventIdxs <- list()
  for (Season in c("S","W")) EventIdxs[[Season]] <- which(External$Season==Season) - 1 # the -1 for C++ indexing
  # exclude the last event as Pcp is trimmed to end of last year - always occurs in winter
  EventIdxs[["W"]] <- EventIdxs[["W"]][1:(length(EventIdxs[["W"]])-1)]

  # Calculate the total volume of small events to be added
  SmallEvents.Volume <- c(S = sum(External$WSA[External$Season=="S"])*Summer.Ratio/100,
                          W = sum(External$WSA[External$Season=="W"])*Winter.Ratio/100)

  # Number of small events to insert
  SmallEvents.N <- c(S = round( unname(SmallEvents.Volume["S"]) / mean(SmallEvents[["S"]]$WSA) ),
                     W = round( unname(SmallEvents.Volume["W"]) / mean(SmallEvents[["W"]]$WSA) ))
  print.noquote(paste("Small events to add (S):", round(SmallEvents.N["S"]/Years), "[N/season],", round((SmallEvents.Volume["S"])/Years,1), "[mm/season]"))
  print.noquote(paste("Small events to add (W):", round(SmallEvents.N["W"]/Years), "[N/season],", round((SmallEvents.Volume["W"])/Years,1), "[mm/season]"))

  # Randomly sample from small events for N
  SmallEvents.Sample <- list()
  for (Season in c("S","W")) {
    # the 1.05 is to deliberately sample more small events than needed - we trim to exact amount next
    SmallEvents.Sample[[Season]] <- SmallEvents[[Season]][ sample(nrow(SmallEvents[[Season]]), size = SmallEvents.N[Season]*1.05, replace = TRUE ) ,  ]
    # Trim SmallEvents.Sample to required volume
    SmallEvents.Sample[[Season]] <- SmallEvents.Sample[[Season]][1:which.max( cumsum(SmallEvents.Sample[[Season]]$WSA) > SmallEvents.Volume[Season]) , ]
    # Convert to time steps
    SmallEvents.Sample[[Season]]$WSD <- SmallEvents.Sample[[Season]]$WSD/60
    # work out WSI, as this is added to each timestep
    SmallEvents.Sample[[Season]]$WSI <- 100*SmallEvents.Sample[[Season]]$WSA / SmallEvents.Sample[[Season]]$WSD # units 1/100th mm
    SmallEvents.Sample[[Season]][,c("WSA","Distance") ] <- NULL # no longer needed
    storage.mode(SmallEvents.Sample[[Season]]$WSD) <- "integer"
    storage.mode(SmallEvents.Sample[[Season]]$WSI) <- "integer"
  }

  print.noquote(paste0("Initial pcp sum: ", round(sum(Pcp)/(Years*100),2), " [mm/yr]" ))

  # Now actually add the small events, performed for each season separately
  for (Season in c("W","S")) {
    Pcp <- insert_small_events_station(Pcp = Pcp, Rain = as.logical(Pcp),
                                      EventIdxs = EventIdxs[[Season]],
                                      EventStart = External$StartDSD[1:(nrow(External)-1)], # exclude last event!
                                      EventL = External$DSD[1:(nrow(External)-1)], # exclude last event!
                                      SmallEvents = as.matrix(SmallEvents.Sample[[Season]]))
  }
  print.noquote(paste0("Final pcp sum: ", round(sum(Pcp)/(Years*100),2), " [mm/yr]" ))

  # Create new data.frame for Internal
  Internal <- data.frame(Pcp = as.integer(Pcp),
                         Year = as.integer(lubridate::year(seq(External$Date[1], by=paste(arm$Aggregation,"mins"), length.out=length(Pcp)))) )

  # Write the data frame to disk
  if (saveObjects) {
    if (!dir.exists("Small_Events/")) dir.create("Small_Events/")
    if (!dir.exists(paste0("Small_Events/",armParTextExt(),"/"))) dir.create(paste0("Small_Events/",armParTextExt(),"/"))
    save(Internal, file=paste0("Small_Events/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"), compress = "xz")
  }

  endTime <- Sys.time()
  #print.noquote(paste("Time taken:", round(endTime-startTime,2), "secs"))

  return(Internal)

}
