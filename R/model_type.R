#' Title
#'
#' @return
#' @export
#'
#' @examples
get_model_type <- function() {
  return( arm$ModelType )
}

#' Title
#'
#' @param ModelType
#'
#' @return
#' @export
#'
#' @examples
set_model_type <- function(ModelType) {
  if (!(ModelType %in% c("Pidoto","Callau"))) stop("ModelType must be one of either 'Pidoto' or 'Callau'")
  arm$ModelType <- ModelType
}

