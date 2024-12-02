#' Fit marginal probabilty distributions
#'
#' Fits marginal probability distributions to rainfall event variables for summer and winter separately.
#'
#' @param StationID a character string of the station to process.
#' @param AddNoise logical, should noise be added to WSD and DSD?
#' @param DateFrom restrict events >= DateFrom. Optional.
#' @param DateTo restrict events <= DateTo. Optional.
#' @param ThresholdNAs numeric, above which proportion of NA values within an event should the event be excluded?
#' @param saveObjects logical, should marginals be saved to disk?
#'
#' @return A list containing the marginal distribution parameters for each event variable for summer and winter.
#' @export
#'
#' @examples
fit_marginals <- function(StationID,
                          AddNoise = TRUE,
                          DateFrom = NULL,
                          DateTo = NULL,
                          ThresholdNAs = 0.1,
                          saveObjects = TRUE,
                          ModelType = get_model_type()
                          ) {

  set_model_type(ModelType)
  # Hard coded parameters
  LowerN <- 50 # Under how many samples should we warn the user of a potentially bad fitting?

  print.noquote(paste("## Fit marginals:", StationID, "##"))

  # load the station event data
  events <- loadRObject(paste0("Events/",armParText(),"/",StationID,".RData"))
  # Crop the time series if desired
  if (!is.null(DateFrom)) events <- events[events$StartDateTime >= DateFrom,] # TODO: check validity of DateFrom
  if (!is.null(DateTo)) events <- events[events$StartDateTime <= DateTo,] # TODO: check validity of DateTo
  # Add noise if desired
  if (AddNoise) {
    events$DSD <- events$DSD + runif(nrow(events), min=-(arm$Aggregation/2), max=(arm$Aggregation/2))
    events$WSD <- events$WSD + runif(nrow(events), min=-(arm$Aggregation/2), max=(arm$Aggregation/2))
  }

  # Create a column which describes the ratio of WSP to WSA
  events$WSP_WSA <- events$WSP/events$WSA
  # Create a column which describes the ratio of WSPT to WSD
  events$WSPT_WSD <- events$WSPT/events$WSD

  returnVals <- list(
    S = list(),
    W = list()
  )

  # Now loop through each CP and perform the parameter estimation
  for (Season in c("S","W")) {

    # subset the data by season
    eventsSub <- events[which(events$Season==Season),]
    print.noquote(paste0(seasonText(Season), " event N: ", nrow(eventsSub)))
    if (nrow(eventsSub) < LowerN) warning(paste0("Low number of events for ",seasonText(Season),"! Fitting performance may be affected."))

    #### Estimate DSD PD ####
    vals <- eventsSub$DSD[eventsSub$DryNAs / (eventsSub$DSD/arm$Aggregation) <= ThresholdNAs ] # filter out events with too many NAs
    DSD_Lmoments <- lmomco::lmoms(vals)
    if (isPidoto()) DSD_PD <- lmomco::parwei(DSD_Lmoments, checklmom = FALSE) # Weibull (3-param)
    if (isCallau()) DSD_PD <- lmomco::parkap(DSD_Lmoments, checklmom = FALSE, snap.tau4 = T) # Kappa (4-param)
    DSD_LBound <- lmomco::qlmomco(0, DSD_PD)
    if (DSD_LBound < arm$DSDmin) warning(paste("DSD lower bound < DSDmin! Values below DSDmin are rounded up.")) # Report if the lower bound is below 0

    #### Estimate WSA PD ####
    vals <- eventsSub$WSA[eventsSub$WetNAs / (eventsSub$WSD/arm$Aggregation) <= ThresholdNAs ] # filter out events with too many NAs
    WSA_Lmoments <- lmomco::lmoms(vals)
    if (isPidoto()) WSA_PD <- lmomco::parwei(WSA_Lmoments, checklmom = FALSE) # Weibull (3-param)
    if (isCallau()) WSA_PD <- lmomco::parkap(WSA_Lmoments, checklmom = FALSE, snap.tau4 = T) # Kappa (4-param)
    WSA_LBound <- lmomco::qlmomco(0, WSA_PD)
    if (WSA_LBound < arm$WSAmin) warning(paste("WSA lower bound < WSAmin! Values below WSAmin are rounded up.")) # Report if the lower bound is below 0

    #### Estimate WSD PD ####
    vals <- eventsSub$WSD[eventsSub$WetNAs / (eventsSub$WSD/arm$Aggregation) <= ThresholdNAs ] # filter out events with too many NAs
    WSD_Lmoments <- lmomco::lmoms(vals)
    WSD_PD <- lmomco::parln3(WSD_Lmoments, checklmom = FALSE) # Log-normal (3-param)
    WSD_LBound <- lmomco::qlmomco(0, WSD_PD)
    # check the lower bound
    if (WSD_LBound < arm$Aggregation) warning(paste("WSD lower bound < Aggregation! Values below Aggregation are rounded up.")) # Report if the lower bound is below 0

    # create return object
    returnVals[[Season]] <- list(
      WSA = WSA_PD$para,
      WSD = WSD_PD$para,
      DSD = DSD_PD$para
    )

    #### Estimate the WSP_WSA ####
    if (isPidoto()) {
      vals <- eventsSub$WSP_WSA[eventsSub$WSD!=arm$Aggregation | # exclude events at the aggregation length
                                eventsSub$WetNAs / (eventsSub$WSD/arm$Aggregation) <= ThresholdNAs  ] # or with too many missing vals
      WSP_WSA_Lmoments <- lmomco::lmoms(vals)
      WSP_WSA_PD <- lmomco::parwei(WSP_WSA_Lmoments, checklmom = FALSE) # Weibull (3-param)
      returnVals[[Season]][["WSP_WSA"]] <- WSP_WSA_PD$para # add to return value
    }

    #### Estimate the WSI ####
    if (isCallau()) {
      vals <- eventsSub$WSI[eventsSub$WetNAs / (eventsSub$WSD/arm$Aggregation) <= ThresholdNAs] # exclude events with too many missing vals
      WSI_Lmoments <- lmomco::lmoms(vals)
      WSI_PD <- lmomco::parkap(WSI_Lmoments, snap.tau4 = T) # Kappa (4-param)
      returnVals[[Season]][["WSI"]] <- WSI_PD$para # add to return value
    }

    #### Estimate the WSP ####
    if (isCallau()) {
      vals <- eventsSub$WSP[eventsSub$WetNAs / (eventsSub$WSD/arm$Aggregation) <= ThresholdNAs] # exclude events with too many missing vals
      WSP_Lmoments <- lmomco::lmoms(vals)
      WSP_PD <- lmomco::pargno(WSP_Lmoments, checklmom = FALSE) # Gen Norm (3-param)
      returnVals[[Season]][["WSP"]] <- WSP_PD$para # add to return value
    }

    #### Save the PD variables to file ####
    if (saveObjects) {
      if (!dir.exists("Marginals/")) dir.create("Marginals/")
      # DSD
      if (!dir.exists("Marginals/DSD/")) dir.create("Marginals/DSD/")
      if (!dir.exists(paste0("Marginals/DSD/",armParTextExt()))) dir.create(paste0("Marginals/DSD/",armParTextExt()))
      save(DSD_PD, file=paste0("Marginals/DSD/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      # WSD
      if (!dir.exists("Marginals/WSD/")) dir.create("Marginals/WSD/")
      if (!dir.exists(paste0("Marginals/WSD/",armParTextExt()))) dir.create(paste0("Marginals/WSD/",armParTextExt()))
      save(WSD_PD, file=paste0("Marginals/WSD/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      # WSA
      if (!dir.exists("Marginals/WSA/")) dir.create("Marginals/WSA/")
      if (!dir.exists(paste0("Marginals/WSA/",armParTextExt()))) dir.create(paste0("Marginals/WSA/",armParTextExt()))
      save(WSA_PD, file=paste0("Marginals/WSA/",armParTextExt(),"/",StationID,"-",Season,".RData"))

      if (isPidoto()) {
        # WSP_WSA
        if (!dir.exists("Marginals/WSP_WSA/")) dir.create("Marginals/WSP_WSA/")
        if (!dir.exists(paste0("Marginals/WSP_WSA/",armParTextExt()))) dir.create(paste0("Marginals/WSP_WSA/",armParTextExt()))
        save(WSP_WSA_PD, file=paste0("Marginals/WSP_WSA/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      }
      if (isCallau()) {
        # WSI
        if (!dir.exists("Marginals/WSI/")) dir.create("Marginals/WSI/")
        if (!dir.exists(paste0("Marginals/WSI/",armParTextExt()))) dir.create(paste0("Marginals/WSI/",armParTextExt()))
        save(WSI_PD, file=paste0("Marginals/WSI/",armParTextExt(),"/",StationID,"-",Season,".RData"))
        # WSP
        if (!dir.exists("Marginals/WSP/")) dir.create("Marginals/WSP/")
        if (!dir.exists(paste0("Marginals/WSP/",armParTextExt()))) dir.create(paste0("Marginals/WSP/",armParTextExt()))
        save(WSP_PD, file=paste0("Marginals/WSP/",armParTextExt(),"/",StationID,"-",Season,".RData"))
      }
    }
  }

  names(returnVals) <- c("Summer","Winter") # rename to friendly names
  return(returnVals)

}
