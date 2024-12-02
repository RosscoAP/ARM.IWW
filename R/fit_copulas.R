#' Title
#'
#' @param StationID
#' @param AddNoise
#' @param DateFrom
#' @param DateTo
#' @param ThresholdNAs
#' @param saveObjects
#' @param ModelType
#'
#' @return
#' @export
#'
#' @examples
fit_copulas <- function(StationID,
                          AddNoise = TRUE,
                          DateFrom = NULL,
                          DateTo = NULL,
                          ThresholdNAs = 0.1,
                          saveObjects = TRUE,
                          ModelType = get_model_type()
                        ) {

  set_model_type(ModelType)
  # Hard coded parameters
  LowerN <- 100 # Under how many samples should we warn the user of a potentially bad fitting?

  print.noquote(paste("## Fit copulas:", StationID, "##"))

  # load the station event data
  events <- loadRObject(paste0("Events/",armParText(),"/",StationID,".RData"))
  # Crop the time series if desired
  if (!is.null(DateFrom)) events <- events[events$StartDateTime >= DateFrom,] # TODO: check validity of DateFrom
  if (!is.null(DateTo)) events <- events[events$StartDateTime <= DateTo,] # TODO: check validity of DateTo
  # Some data prep.
  events$WSD_orig <- events$WSD # save a copy of the unmodified WSD, as it will be needed later on
  # Add noise if desired
  if (AddNoise) events$WSD <- events$WSD + runif(nrow(events), min=-(arm$Aggregation/2), max=(arm$Aggregation/2))

  # Do not allow WSP to be greater than WSA
  events$WSP[which(events$WSP>events$WSA)] <- events$WSA[which(events$WSP>events$WSA)]
  # Recalculate the WSI based on the noisy versions
  events$WSI <- events$WSA/events$WSD
  # Create a column representing the ration between WSP and WSA
  events$WSP_WSA <- events$WSP/events$WSA
  events$WSPT_WSD <- events$WSPT/events$WSD

  # create return object
  returnVals <- list(
    S = list(),
    W = list()
  )
  if (saveObjects) {
    if (!dir.exists("Copulas/")) dir.create("Copulas/")
    if (!dir.exists("Copulas/WSI/")) dir.create("Copulas/WSI/")
    if (!dir.exists("Copulas/WSP/")) dir.create("Copulas/WSP/")
  }

  # Now loop through each CP and perform the parameter estimation
  for (Season in c("S","W")) {

    # subset the data by season
    eventsSub <- events[which(events$Season==Season),]
    print.noquote(paste0(seasonText(Season), " event N: ", nrow(eventsSub)))
    if (nrow(eventsSub) < LowerN) warning(paste0("Low number of events for ",seasonText(Season),"! Fitting performance may be affected."))

    if (isPidoto()) {

      #### Estimate WSI (WSA-WSD) Copula ####
      pseudoObs_WSA_WSD <- copula::pobs(eventsSub[,c("WSA","WSD")]) # Create pseudo obs
      copulaWSI <- copula::fitCopula(copula::khoudrajiCopula(copula2 = copula::gumbelCopula(),
                                      shapes = copula::fixParam(c(NA_real_, 1), c(FALSE, TRUE))), # second shape parameter is fixed to 1
                                      data = pseudoObs_WSA_WSD, start = c(1.6, 0.9), optim.method = "Nelder-Mead")@copula
      returnVals[[Season]][["WSI"]] <- c(Alpha = copulaWSI@copula2@parameters, Shape = copulaWSI@shapes[1])

      #### Estimate WSP (WSD-WSP:WSA) Copula ####
      # Create pseudo obs
      pseudoObs_WSD_WSP_WSA <- copula::pobs(eventsSub[(eventsSub$WSD_orig!=arm$Aggregation),c("WSD","WSP_WSA")], ties.method = "random") # excludes the single timestep events!
      pseudoObs_WSD_WSP_WSA[,"WSD"] <- 1 - pseudoObs_WSD_WSP_WSA[,"WSD"] # flip the Y axis so we are in the standard Copula rotation
      copulaWSP <- copula::fitCopula(copula::khoudrajiCopula(copula2 = copula::normalCopula(),
                                     shapes = copula::fixParam(c(1, NA_real_), c(TRUE, FALSE))), # first shape parameter is fixed to 1
                                     data = pseudoObs_WSD_WSP_WSA, start = c(0.7, 0.9), optim.method = "Nelder-Mead")@copula
      returnVals[[Season]][["WSP"]] <- c(Rho = copulaWSP@copula2@parameters, Shape = copulaWSP@shapes[2])

      #### Save the copula objects to file ####
      if (saveObjects) {
        # WSI
        if (!dir.exists(paste0("Copulas/WSI/",armParTextExt()))) dir.create(paste0("Copulas/WSI/",armParTextExt()))
        save(copulaWSI, file=paste0("Copulas/WSI/",armParTextExt(),"/",StationID,"-",Season,".RData"))
        # WSP
        if (!dir.exists(paste0("Copulas/WSP/",armParTextExt()))) dir.create(paste0("Copulas/WSP/",armParTextExt()))
        save(copulaWSP, file=paste0("Copulas/WSP/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      }

    }

    if (isCallau()) {

      #### Estimate WSP (WSI:WSP) Copula ####
      pseudoObs_WSI_WSP <- copula::pobs(eventsSub[,c("WSI","WSP")]) # Create pseudo obs
      copulaWSP <- copula::fitCopula(copula::normalCopula(), data = pseudoObs_WSI_WSP)@copula
      returnVals[[Season]][["WSP"]] <- c(Rho = copulaWSP@parameters)

      #### Estimate WSI (WSD:WSA) Copula ####
      # In the regional empirical Copula, the WSA and WSD are normalised with their median
      # Conversion into the unit space occurs during synthesis!
      copulaWSI <- cbind(WSA = eventsSub$WSA / median(eventsSub$WSA), WSD = eventsSub$WSD / median(eventsSub$WSD))

      #### Save the copula objects to file ####
      if (saveObjects) {
        # WSI
        if (!dir.exists(paste0("Copulas/WSI/",armParTextExt()))) dir.create(paste0("Copulas/WSI/",armParTextExt()))
        save(copulaWSI, file=paste0("Copulas/WSI/",armParTextExt(),"/",StationID,"-",Season,".RData"))
        # WSP
        if (!dir.exists(paste0("Copulas/WSP/",armParTextExt()))) dir.create(paste0("Copulas/WSP/",armParTextExt()))
        save(copulaWSP, file=paste0("Copulas/WSP/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      }

    }
  }

  names(returnVals) <- c("Summer","Winter") # rename to friendly names
  return(returnVals)

}
