#' Title
#'
#' @param StationID
#' @param Years
#' @param StartYear
#' @param Realisation
#' @param RegEmpCopula
#' @param saveObjects
#' @param ModelType
#'
#' @return
#' @export
#'
#' @examples
gen_external <- function(StationID,
                        Years,
                        StartYear = 1900,
                        Realisation = NA,
                        RegEmpCopula = StationID,
                        saveObjects = TRUE,
                        ModelType = get_model_type()
                        ) {

  # TODO: check parameters
  set_model_type(ModelType)

  # Hard coded parameters
  # Roughly how many samples should be generated in order to fill the simulation period?
  numSamples <- 400*Years # better higher than lower

  print.noquote(paste("## Generate external time series: ", StationID, "##"))
  startTime <- Sys.time()

  # work out start and end dates
  StartDate <- as.POSIXct(paste0(StartYear,"-01-01"), tz="UTC")
  EndDate   <- as.POSIXct(paste0(StartYear+Years,"-01-01"), tz="UTC")
  # For efficiency, use timesteps instead of dates during the synthesis, and only convert to dates at the end
  totalTimeSteps <- as.integer((as.numeric(EndDate)-as.numeric(StartDate))/(60*arm$Aggregation))

  # List to hold the WSA, WSD, WSD samples
  WSA <- list()
  WSD <- list()
  DSD <- list()
  WSD_Pvals <- list()
  WSA_WSD <- list() # holds the samples in the [0,1] space

  for (Season in c("S","W")) {

    # Load the marginal distributions
    PD_DSD <- loadRObject(file=paste0("Marginals/DSD/",armParTextExt(),"/",StationID,"-",Season,".RData"))
    PD_WSA <- loadRObject(file=paste0("Marginals/WSA/",armParTextExt(),"/",StationID,"-",Season,".RData"))
    PD_WSD <- loadRObject(file=paste0("Marginals/WSD/",armParTextExt(),"/",StationID,"-",Season,".RData"))

    # Load Copulas
    if (isPidoto()) { # Khoudraji Gumbel Copula
      copulaWSI <- loadRObject(file=paste0("Copulas/WSI/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      # Generate n number of samples from the Copula
      WSA_WSD[[Season]] <- copula::rCopula(ceiling(numSamples*0.5), copulaWSI) # 0.5 as each season is half the year
    }
    if (isCallau()) { # Regional Empirical Copula
      # Matrix to hold copula values
      WSA_WSD[[Season]] <- matrix(data = NA, nrow = 0, ncol = 2)
      # Load empirical Copulas from each station from disk
      for (i in unique(c(StationID,RegEmpCopula)) ) WSA_WSD[[Season]] <- rbind(WSA_WSD[[Season]], loadRObject(paste0("Copulas/WSI/",armParTextExt(),"/",StationID,"-",Season,".RData")) )
      # Note for previous line: unique(c(StationID,RegEmpCopula)) ensures StationID is included but not duplicated if included in RegEmpCopula
      # Generate pseudo obs (as until now they are just normalised values)
      WSA_WSD[[Season]] <- copula::pobs(WSA_WSD[[Season]])
      # and randomly sample with replacement
      WSA_WSD[[Season]] <- WSA_WSD[[Season]][ sample(x = nrow(WSA_WSD[[Season]]), size = ceiling(numSamples*0.5), replace = T),  ] # 0.5 as each season is half the year
    }

    # And then transform from the unit space using the marginals
    WSD_Pvals[[Season]] <- WSA_WSD[[Season]][,2] # we need to save the pvals as they are used in script 4,
    # in the past we used the plmomco() in script 4, but as this was applied after rounding it caused issues where WSDs were below the lower bound of the PD)
    WSD[[Season]] <- lmomco::qlmomco(WSA_WSD[[Season]][,2], PD_WSD)
    WSA[[Season]] <- lmomco::qlmomco(WSA_WSD[[Season]][,1], PD_WSA)
    # Generate also DSD samples which are independent from wetspells
    DSD[[Season]] <- lmomco::rlmomco(ceiling(numSamples*0.5), PD_DSD) # 0.5 as each season is half the year
    rm(copulaWSI) # object no longer needed

    # Round the times to the specified output timestep
    WSD[[Season]] <- round(WSD[[Season]]/arm$Aggregation,0)*arm$Aggregation
    DSD[[Season]] <- round(DSD[[Season]]/arm$Aggregation,0)*arm$Aggregation

    # Make sure that the values are above the required lower bounds - this can occur due to problems with PD fittings
    if (min(WSD[[Season]])<0) warning(paste("Negative WSD found! Values rounded to timestep length."))
    if (min(WSA[[Season]])<arm$WSAmin) warning(paste("WSA<WSAmin found! Values rounded to WSAmin."))
    DSD[[Season]][which(DSD[[Season]]<arm$DSDmin)] <- arm$DSDmin
    WSD[[Season]][which(WSD[[Season]]<arm$Aggregation)] <- arm$Aggregation
    WSA[[Season]][which(WSA[[Season]]<arm$WSAmin)] <- arm$WSAmin

  }
  rm(WSA_WSD) # object no longer needed
  #print.noquote(paste(Sys.time(),"- Copula/marginal samples generated, starting event generation."))

  # Create blank vectors to store the relevant time series data
  # The vectors will be the length of all possible events, and then later trimmed to the number of actual events
  # To be honest I don't know why I just don't use directly a data.frame, maybe there was a good reason to do it like this at the time.
  TS.startDateTime <- rep(NA,numSamples) # start datetime of the event
  TS.year <- rep(NA,numSamples) # YYYY
  TS.wsa <- rep(NA,numSamples)  # [mm]
  TS.wsd <- rep(NA,numSamples)  # [min]
  TS.wsd_pval <- rep(NA,numSamples)  # [0-1]
  TS.dsd <- rep(NA,numSamples)  # [min]
  TS.wsi <- rep(NA,numSamples)  # [mm/min]
  TS.Season_dry <- rep(NA,numSamples) # Season number for dry
  TS.Season_wet <- rep(NA,numSamples) # Season number for wet

  # Reset the current timestep
  currentTimeStep <- 0L
  currentDate <- StartDate
  idx <- 0L # initialise idx counter
  iWet <- iDry <- c(S = 1, W = 1)  # Keep counters of the next to be sampled wet/dry event

  # Loop until we are > than EndDate
  while (currentTimeStep < totalTimeSteps) {
    idx <- idx + 1L  # increase current event index
    for (WetDry in c("Wet","Dry")) { # Wet spell and dry spell chosen alternating
      iSeason <- returnSeason(currentDate) # season the current event falls in
      # Now add the randomly chosen event to the time series
      if (WetDry=="Wet") {
        TS.wsd[idx] <- WSD[[iSeason]][ iWet[iSeason] ]
        TS.wsd_pval[idx] <- WSD_Pvals[[iSeason]][ iWet[iSeason] ]
        TS.wsa[idx] <- WSA[[iSeason]][ iWet[iSeason] ]
        currentTimeStep <- currentTimeStep + as.integer(TS.wsd[idx]/arm$Aggregation)
        currentDate <- currentDate + as.integer(TS.wsd[idx])*60
        iWet[iSeason] <- iWet[iSeason] + 1L  # increase count by 1
      } else {
        TS.dsd[idx] <- DSD[[iSeason]][ iDry[iSeason] ]
        currentTimeStep <- currentTimeStep + as.integer(TS.dsd[idx]/arm$Aggregation)
        currentDate <- currentDate + as.integer(TS.dsd[idx])*60
        iDry[iSeason] <- iDry[iSeason] + 1L  # increase count by 1
      }
    }
  }

  # Calculate the WSI
  TS.wsi<- round(TS.wsa/TS.wsd,5)
  TS.wsa<- round(TS.wsa,5) # round the WSA
  TS.startDateTime <- StartDate + ((cumsum(TS.wsd) + cumsum(TS.dsd))*60) - ((TS.wsd + TS.dsd)*60)
  TS.year <- lubridate::year(TS.startDateTime)

  # Create a data frame for export
  External <- data.frame(
    Date = as.POSIXct(TS.startDateTime[1:idx], origin="1970-01-01", tz="UTC"),
    Year = TS.year[1:idx],
    Season = returnSeason(as.POSIXct(TS.startDateTime[1:idx], origin="1970-01-01", tz="UTC")),
    WSD = TS.wsd[1:idx],
    WSD_Pval = TS.wsd_pval[1:idx],
    DSD = TS.dsd[1:idx],
    WSA = TS.wsa[1:idx],
    WSI = TS.wsi[1:idx]
  )

  # Write the data frame to disk
  if (saveObjects) {
    if (!dir.exists("External/")) dir.create("External/")
    if (!dir.exists(paste0("External/",armParTextExt(),"/"))) dir.create(paste0("External/",armParTextExt(),"/"))
    save(External, file=paste0("External/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"))
  }

  endTime <- Sys.time()
  print.noquote(paste("Time taken:", round(endTime-startTime,2), "secs"))

  return(External)

}
