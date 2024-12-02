#' Title
#'
#' @return
#' @export
#'
#' @examples
get_summer_months <- function() {
  return( arm$SummerMonths )
}


#' Title
#'
#' @param Months
#'
#' @return
#' @export
#'
#' @examples
set_summer_months <- function(Months) {
  # need to include test that Months only contains integers in range 1-12
  arm$SummerMonths <- Months
}

