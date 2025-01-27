#' Title
#'
#' @param StationID
#' @param Years
#' @param Realisation
#' @param Alpha
#' @param saveSpatialObjects
#' @param saveObjects
#' @param ModelType
#'
#' @return
#' @export
#'
#' @examples
gen_internal <- function(StationID,
                        Years,
                        Realisation = NA,
                        Alpha = 1/3,
                        saveSpatialObjects = TRUE,
                        saveObjects = TRUE,
                        ModelType = get_model_type()
                        ) {

  set_model_type(ModelType)

  # Hard coded parameter
  # At what resolution will the internal structure be produced before being aggregated to the output time step?
  # A lower number will lead to lower rounding errors
  InternalResolution <- 1 # [min]
  nSamples <- 1e6 # how many samples should be generated for the numerical implementation of cCopula

  print.noquote(paste("## Generate internal time series:", StationID, "##"))
  startTime <- Sys.time()

  # Load the external structure of the time series
  Events <- loadRObject(paste0("External/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"))
  totalLength <- (sum(Events$DSD)+sum(Events$WSD))/arm$Aggregation # in timesteps

  # Create blank columns within Events
  Events$WSP <- NA
  Events$WSPT <- NA
  Events$Xi <- NA
  Events$Pcp <- as.list(NA)
  Events$Season <- returnSeason(Events$Date)

  if (isCallau()) { # Callau version of the model

    # this column is only needed in the Callau model
    Events$WSI_CummProb <- 0

    for (Season in c("S","W")) { # loop around seasons
      # Load objects
      WSP_Copula <- loadRObject(paste0("Copulas/WSP/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      WSP_PD <- loadRObject(paste0("Marginals/WSP/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      WSI_PD <- loadRObject(paste0("Marginals/WSI/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      # In the Callau model, we calculate back from WSI
      WSI_CummProb <- lmomco::plmomco(Events$WSI[Events$Season==Season & Events$WSD!=arm$Aggregation], WSI_PD)
      pVals <- copula::cCopula(cbind((WSI_CummProb), runif(length(WSI_CummProb))), copula = WSP_Copula, inverse = T, indices = 2)
      Events$WSP[Events$Season==Season & Events$WSD!=arm$Aggregation] <- mapply(FUN = max, lmomco::qlmomco(pVals, WSP_PD), 1e-4) # Do not allow a negative value to be returned.
    }
    # and events at the aggregation level
    Events$WSP[Events$WSD==arm$Aggregation] <- Events$WSA[Events$WSD==arm$Aggregation]
    # Due to the WSP calculation process, it can be that the volume of the WSP time step is greater than WSA.
    # As such, ensure that WSP<=WSA
    Events$WSP[which(Events$WSP>Events$WSA)] <- Events$WSA[which(Events$WSP>Events$WSA)]
  }
  if (isPidoto()) { # Pidoto version of the model

    # The Pidoto model uses a temporary workaround method for Khoudraji type Copulas
    # As copula::cCopula method isn't implemented yet. So we use a numerical method generated first random samples,
    # and then the closest sample to the required pValue is returned.
    # See https://stackoverflow.com/questions/43472234/r-fastest-way-to-find-nearest-value-in-vector

    # this column is only needed in the Pidoto model
    Events$WSP_WSA <- NA

    # Loop around Seasons (and seasons) and generate a large number of sample Copula values
    for (Season in c("S","W")) {
      # Load the WSP/WSI Copula and marginal distributions
      WSP_Copula <- loadRObject(paste0("Copulas/WSP/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      WSP_WSA_PD <- loadRObject(paste0("Marginals/WSP_WSA/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      # generate randome samples across the entire 2D unit space
      pSim <- data.table::as.data.table(copula::rCopula(nSamples, WSP_Copula))
      colnames(pSim) <- c("U","V") # U: WSD, V; WSP_WSA
      data.table::setkeyv(pSim,"U")
      IdxVals <- which(Events$WSD!=arm$Aggregation & Events$Season==Season) # all events in current season > Aggregation
      pValsWSD <- data.table::data.table(U = (1-Events$WSD_Pval[IdxVals]), Idx = IdxVals ) # we also need to store the idx values becuase when the merge is performed it will automatically sort the data.table
      data.table::setkeyv(pValsWSD,"U") # create a data.table index (key)
      pValsWSP <- pSim[pValsWSD, roll='nearest'] # relies on .datatable.aware being true
      # Then transform using the marginal distribution
      Events$WSP_WSA[pValsWSP$Idx] <- lmomco::qlmomco(pValsWSP$V, WSP_WSA_PD)
    }
    rm(pValsWSD, pValsWSP, pSim, IdxVals) # objects no longer needed
    gc()
    Events$WSP_WSA[which(Events$WSP_WSA < 0.05)] <- 0.05 # Do not allow too small values (obs generally always > 0.05)
    Events$WSP_WSA[which(Events$WSP_WSA > 0.995)] <- 0.995 # Nor a value above 1
    Events$WSP_WSA[(Events$WSD==arm$Aggregation)] <- 1 # Force the value to be 1 if WSD==Aggregation
    Events$WSP <- Events$WSP_WSA * Events$WSA # Then multiply by WSA to get the WSP
  }

  Events$WSP <- Events$WSP/arm$Aggregation # To match up with the other methods, the WSP is expressed as an intensity
  # In order to preserve our WSP (over a full output time step), we will manually insert the WSP after calculating the exponential function
  # This has the following ramifications:
  # - for our exponential function, WSD = WSD - arm$Aggregation
  # - the WSPT is relative to this modified WSD
  # - for our exponential function, WSA = WSA - (WSP*arm$Aggregation)
  Events$WSD_mod <- Events$WSD - arm$Aggregation # modified WSD for use in the exponential function
  Events$WSA_mod<- abs(round(Events$WSA-(Events$WSP*arm$Aggregation),5)) # modified WSA for use in the exponential function
  # Determine the WSPT (wet spell peak time) for each event
  # According to Haberlandt (1998), value is chosen from a uniform distribution
  Events$WSPT <- round(  sapply(X=Events$WSD_mod, FUN=runif, min=0, n=1)/arm$Aggregation,0)*arm$Aggregation

  # Functions to first find Xi and then also return parts of the P equation
  fnXi <- function(x, WSP, WSPT, WSD, WSA) return (((WSP/x) * (2 - exp(-x*WSPT) - exp(x*(WSPT-WSD) ))) - WSA)
  findXi <- function(WSP, WSPT, WSD, WSA)  return (uniroot(fnXi, WSP=WSP, WSPT=WSPT, WSD=WSD, WSA=WSA, lower=1e-4, upper=10, tol=1e-5, extendInt="yes", maxiter=1e3)$root)
  cPart <- function(WSPT,WSD) return (c( rep(1,round(WSPT/InternalResolution,0)), rep(-1,round((WSD-WSPT)/InternalResolution,0))) ) # This is the c parameter used within the Pcp formula (see Haberlandt (1998))
  tPart <- function(WSPT,WSD) return ( seq(from=(InternalResolution-WSPT), to=(WSD-WSPT), by=InternalResolution) ) # the time part of the function
  # Function to calculate pcp intensity at each timestep
  returnPcp <- function(WSP,Xi,WSPT,WSD) {
    if (WSD==0) return(  ) # return nothing if duration is zero
    if (is.na(Xi)) return( rep(0,times=WSD/arm$Aggregation) ) # return a series of zeroes if Xi is NA, which means there is no exponential part
    return( # otherwise return the double exponential function
      unname(round(
        tapply( # use tapply() to aggregate up to the output time step
          X= (WSP * exp( cPart(WSPT,WSD) * Xi * tPart(WSPT,WSD) )),
          INDEX=rep(c(1:(WSD/arm$Aggregation)),each=(arm$Aggregation/InternalResolution)), FUN=mean)*arm$Aggregation
        ,4)
      ))
  }
  # Function which inserts the peak volume in the right place in order to preserve the WSP volume
  insertPeak <- function(Pcp,WSPT,WSP) {
    if (WSPT==0) return( c(round(WSP*arm$Aggregation, 4), Pcp)  ) # Peak occurs at the start
    if (WSPT==(length(Pcp)*arm$Aggregation)) return( c(Pcp, round(WSP*arm$Aggregation, 4))  ) # Peak occurs at the end
    # Otherwise peak occurs in the middle
    return( c(Pcp[1:(WSPT/arm$Aggregation)], round(WSP*arm$Aggregation, 4), Pcp[((WSPT/arm$Aggregation)+1):length(Pcp)]) )
  }

  # these are the events which need to be solved for Xi (i.e. (WSD-arm$Aggregation) > 0 AND WSA > 0)
  toBeSolved <- which(Events$WSD_mod > 0 & Events$WSA_mod > 1e-4)
  Events$WSI_mod <- Events$WSI #(Events$WSA-(Events$WSP*arm$Aggregation))/(Events$WSD-arm$Aggregation) # modified WSI considering one less timestep and less volume
  Events$WSP_mod <- Events$WSI_mod + (Alpha * (Events$WSP - Events$WSI_mod)  )
  Events$Xi[toBeSolved] <- mapply(FUN=findXi, WSP=Events$WSP_mod[toBeSolved], WSPT=Events$WSPT[toBeSolved], WSD=Events$WSD_mod[toBeSolved], WSA=Events$WSA_mod[toBeSolved])
  Events$Pcp <- mapply(FUN=returnPcp, WSPT=Events$WSPT, WSD=Events$WSD_mod, Xi=Events$Xi, WSP=Events$WSP_mod)

  # Insert the preserved WSP at the right location
  Events$Pcp <- mapply(FUN=insertPeak,Pcp=Events$Pcp,WSPT=Events$WSPT,WSP=Events$WSP)

  Events$SumPcp <- sapply(FUN=sum,Events$Pcp) # for debugging purposes
  Events$Error <- Events$SumPcp/Events$WSA # error per event
  OverallWSAError <- round( (sum(Events$SumPcp)-sum(Events$WSA))*100/sum(Events$WSA) ,2)
  #print.noquote(paste("Overall WSA error (%): ",OverallWSAError)) # for debugging

  # Write the data frame to disk
  if (saveObjects) {
    # Save the external structure again to disk with extra information needed for wet spell peaks
    External <- Events[,c("Date","Year","Season","WSD","WSD_Pval","DSD","WSA","WSI","WSP","WSPT")] # to conform with previous naming
    save(External, file=paste0("External/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"))
    rm(External) # remove temporary object
    ##### Create data.frame to store output time series
    Internal <- data.frame(
      Pcp = double(totalLength), # Set Pcp first to zero. This will account for the DSDs. Later the WSDs will be overwritten
      Year = as.integer(lubridate::year(seq(Events$Date[1], by=paste(arm$Aggregation,"mins"), length.out=totalLength)))
    )
    # Determine the start and end indexes for each wet event
    Events$startIndex <- (as.numeric(Events$Date)-as.numeric(Events$Date[1]))/(60*arm$Aggregation)+1
    Events$endIndex <- Events$startIndex + (Events$WSD/arm$Aggregation) - 1
    # Copy Pcp events from events table to Pcp ts output
    Internal$Pcp[unlist(mapply(FUN=":", Events$startIndex, Events$endIndex))] <- unlist(Events$Pcp)
    # remove any timesteps which extend beyond the year range
    Internal <- Internal[(Internal$Year <= lubridate::year(Events$Date[1])+Years),]
    # Convert Pcp to 1/100th of a mm
    Internal$Pcp <- as.integer(Internal$Pcp*100)
    if (!dir.exists("Internal/")) dir.create("Internal/")
    if (!dir.exists(paste0("Internal/",armParTextExt(),"/"))) dir.create(paste0("Internal/",armParTextExt(),"/"))
    save(Internal, file=paste0("Internal/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"), compress = "xz")
  }

  # Save additionally required output objects for extension into space
  if (saveSpatialObjects) {
    # work out season for dry and wet spells separately
    Events$Season_Wet <- Events$Season
    Events$Season_Dry <- returnSeason(Events$Date + Events$DSD*60)
    # Create a copy of Events data.frame with only the columns needed
    EventsSA <- Events[, c("Season_Wet","Season_Dry","Year","Pcp","DSD","WSA","WSD","WSI","WSP","WSPT","Date")]
    EventsSA$PcpDSD <- sapply(EventsSA$DSD/60,FUN=rep,x=0) # create vectors of zeros for the DSD part
    colnames(EventsSA)[match("Pcp",colnames(EventsSA))] <- "PcpWSD" # rename the Pcp col
    EventsSA$EventIndex <- c(1:nrow(EventsSA)) # simply the row number, so that shuffled timeseries can be recreated
    if (!dir.exists("Internal_Space/")) dir.create("Internal_Space/")
    if (!dir.exists(paste0("Internal_Space/",armParTextExt(),"/"))) dir.create(paste0("Internal_Space/",armParTextExt(),"/"))
    save(EventsSA, file=paste0("Internal_Space/",armParTextExt(),"/",StationID,"-",ifelse(is.na(Realisation),"",paste0(pad0s(Realisation),"-")),Years,"yrs.RData"), compress = "xz")
  }

  endTime <- Sys.time()
  print.noquote(paste("Time taken:", round(endTime-startTime,2), "secs"))

  return(Internal)

}
