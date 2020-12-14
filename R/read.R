#'
#'
#' @export
readMaf <- function(x, ecol = NA, ...) {

  data <- read.table(x, sep = "\t", header = TRUE)

  # rename column names according to conventions
  names(data)[names(data) == "mass_to_charge"] <- "mz"
  names(data)[names(data) == "retention_time"] <- "rtime"

  # check and convert rtime
  if(!any(data$rtime > 60)) {

    data$rtime <- data$rtime * 60

  }

  # add ID if not present
  if(!"id" %in% names(data)) {

    data <- cbind(id = .createId(data), data)

  }

  if(length(ecol == 1) && is.na(ecol)) {
    ecol <- 23:ncol(data)
  }

  # create QFeatures
  readQFeatures(data, ecol = ecol, ...)

}

#'
#'
#' place holder for mzTab
readmzTab <- function(x) {

}

.createId <- function(x) {

  if(!all(c("mz", "rtime") %in% colnames(x))) {

    stop("missing mz and rtime information")

  }

  if("ccs" %in% colnames(x)) {

    idlength <- nchar(nrow(x))

    ids <- paste0("FT",
                  sprintf(eval(paste0("%0", idlength, "d")), row(x)),
                  "M",
                  x$mz,
                  "T",
                  as.integer(x$rtime), # add formating on rtime
                  "C",
                  x$ccs)

  } else {

    idlength <- nchar(nrow(x))

    ids <- paste0("FT",
                  sprintf(eval(paste0("%0", idlength, "d")), row(x)),
                  "M",
                  x$mz,
                  "T",
                  x$rtime)


  }

  return(ids)

}
