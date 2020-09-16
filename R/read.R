#'
#'
#' @export
readMaf <- function(x, ecol = NA, ...) {

  data <- read.table(x, sep = "\t", header = TRUE)

  # rename column names according to conventions
  names(data)[names(data) == "mass_to_charge"] <- "mz"
  names(data)[names(data) == "retention_time"] <- "rtime"

  # check and convert rtime
  # TODO

  if(length(ecol == 1) && is.na(ecol)) {
    ecol <- 22:ncol(data)
  }

  # create QFeatures
  readQFeatures(data, ecol = ecol, ...)


}
