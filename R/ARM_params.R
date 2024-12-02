#' Title
#'
#' @return
#' @export
#'
#' @examples
get_ARM_pars <- function() {
  return( list(Aggregation = arm$Aggregation, DSDmin = arm$DSDmin, WSAmin = arm$WSAmin) )
}

#' Title
#'
#' @param Aggregation
#' @param WSAmin
#' @param DSDmin
#'
#' @return
#' @export
#'
#' @examples
set_ARM_pars <- function(Aggregation = NULL, WSAmin = NULL, DSDmin = NULL) {
  if (!is.null(Aggregation)) arm$Aggregation <- Aggregation
  if (!is.null(WSAmin)) arm$WSAmin <- WSAmin
  if (!is.null(DSDmin)) arm$DSDmin <- DSDmin
}
