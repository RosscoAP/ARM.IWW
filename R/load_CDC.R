#' Title
#'
#' @param StationID
#' @param sourceDir
#' @param includeRecent
#'
#' @return
#' @export
#'
#' @examples
load_CDC <- function(StationID, sourceDir, includeRecent = FALSE) {

  # retreive list of all files in source directory
  sourceFiles <- dir(paste0(sourceDir,"historical/"), pattern = ".zip")
  # Which file within sourceFiles is for this station
  StationFilename <- sourceFiles[ which(grepl(paste0("stundenwerte_RR_",pad0s(as.numeric(StationID),5)), sourceFiles)) ]
  if (length(StationFilename)>1) stop(paste("More than one input file for station",StationID,"found."))
  if (length(StationFilename)!=1) stop(paste("Input file for station",StationID,"not found."))
  # split the filename into parts so we can construct the filename within the ZIP file
  temp <- unlist(strsplit(StationFilename[1],"_")) # 3: station id, 4: date from, 5: date to
  #print.noquote(paste("Station",StationID,"from: ", as.Date(temp[4], format = "%Y%m%d" ), "to: ", as.Date(temp[5], format = "%Y%m%d" )))
  # Unzip the file
  ExtractedFile <- unzip(paste0(sourceDir,"historical/", StationFilename[1]), files = paste0("produkt_rr_stunde_",temp[4],"_",temp[5],"_",temp[3],".txt") )
  Rainfall <- read.csv(ExtractedFile, sep = ";")
  file.remove(ExtractedFile) # delete temporary file
  # also include recent data is requested
  if (includeRecent) { # process is slightly different to above, as different naming scheme
    # only do if recent file exists
    if (file.exists(paste0(sourceDir,"recent/", "stundenwerte_RR_",pad0s(as.numeric(StationID), 5),"_akt.zip"))) {
      # first get the filename we want to extract
      ZipFiles <- unzip(paste0(sourceDir,"recent/", "stundenwerte_RR_",pad0s(as.numeric(StationID), 5),"_akt.zip"), list = T)$Name
      # then extract
      ExtractedFileRecent <- unzip(paste0(sourceDir,"recent/", "stundenwerte_RR_",pad0s(as.numeric(StationID), 5),"_akt.zip"),
                                   files = ZipFiles[ grep("produkt_rr_stunde", ZipFiles) ] )
      RainfallRecent <- read.csv(ExtractedFileRecent, sep = ";")
      file.remove(ExtractedFileRecent) # delete temporary file
      # and merge with Rainfall data.frame - but only dates that aren't already in there
      Idxs <- !(RainfallRecent$MESS_DATUM %in% Rainfall$MESS_DATUM)
      Rainfall <- rbind(Rainfall, RainfallRecent[Idxs,])
    }
  }
  # clean up
  Rainfall$MESS_DATUM <- as.character(Rainfall$MESS_DATUM)
  Rainfall$MESS_DATUM <- as.POSIXct(Rainfall$MESS_DATUM, tz = "UTC", tryFormats = c("%Y%m%d%H"))
  # Make sure the time series starts at 12:00
  if (match(0, as.POSIXlt(Rainfall$MESS_DATUM)$hour)!=1) Rainfall <- Rainfall[match(0, as.POSIXlt(Rainfall$MESS_DATUM)$hour):nrow(Rainfall),]
  #Rainfall$Diff <- c(0,diff(as.numeric(Rainfall$MESS_DATUM))/3600) # difference in hours between time steps
  Rainfall$RS_IND <- as.logical(Rainfall$RS_IND)
  # Create a column which specifies missing values
  Rainfall$isNA <- F; Rainfall$isNA[Rainfall$R1==-999] <- T
  Rainfall$R1[Rainfall$R1==-999] <- 0 # remove -999 values
  # remove unneeded columms
  Rainfall <- Rainfall[,c("MESS_DATUM", "R1", "RS_IND", "isNA")]
  # make sure the time series is continuous - no gaps between hourly values.
  # To do this we will do a join with a blank time series
  Rainfall <- merge(x = Rainfall, y = data.frame(MESS_DATUM = seq.POSIXt(from = min(Rainfall$MESS_DATUM), to = max(Rainfall$MESS_DATUM), by = "hour" )),
                    all.y = T, by = c("MESS_DATUM"))
  Rainfall$isNA[is.na(Rainfall$R1)] <- T
  Rainfall$R1[is.na(Rainfall$R1)] <- 0
  Rainfall$RS_IND[is.na(Rainfall$RS_IND)] <- F
  colnames(Rainfall) <- c("DateTime","Pcp","Rain","isNA")
  #if (!honourRS_IND) Rainfall$Rain <- as.logical(Rainfall$Pcp) # if we don't honour RS_IND, Rain indictor is TRUE where Pcp > 0
  print.noquote(paste0("Loaded station ",StationID,
                       ". Date range: ", min(Rainfall$DateTime), " to ", max(Rainfall$DateTime),
                       ". Missing vals: ", round(100*sum(Rainfall$isNA)/nrow(Rainfall), 2),"%" ))

  return(Rainfall)

}
