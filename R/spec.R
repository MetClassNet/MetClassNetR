#'
#'
#' @export
fillSpectra <- function(QFeatures, Spectra) {

  checkQFeatures(QFeatures)
  checkSpectra(Spectra)

  # switch backends
  Spectra <- setBackend(Spectra, MsBackendMemory())

  ms1_data <- rowData(QFeatures[[1]])
  ms2_data <- spectraData(Spectra, c("id"))

  diff_ids <- setdiff(ms1_data$id, ms2_data$id)

  for(id in diff_ids) {

    # get information
    precursorMz <- ms1_data[ms1_data$id == id,"mz"]
    rtime <- ms1_data[ms1_data$id == id,"rtime"]

    spd <- DataFrame(
      msLevel = c(2L),
      polarity = c(NA_integer_),
      precursorMz = as.numeric(precursorMz),
      rtime = rtime,
      id = id)

    ## Assign m/z and intensity values.
    spd$mz <- list(c())
    spd$intensity <- list(c())

    sps <- Spectra::Spectra(spd, backend = MsBackendDataFrame())

    Spectra <- c(Spectra, sps)

  }

  Spectra

}


#' @title Calculate molecular network from MS2 data
#'
#' @name spec_molNetwork
#'
#' @description
#'
#' `spec_molNetwork` uses `Spectra` objects and calculates molecular similarity
#'     networks based on the comparsion method(s) defined in `methods`
#'
#' @param x `Spectra` Spectra object with single spectrum per feature
#' @param methods `character` Name of function(s) used for comparison
#' @param ... Parameters used by compareSpectra and sub functions
#'
#' @return `list` a named list with adjency matrices for each method
#'
#' @author Michael Witting
#'
#' @export
#'
spec_molNetwork <- function(x, methods = c("ndotproduct"), ...) {

  # sanity checks
  if(!"id" %in% spectraVariables(x)) {

    stop("Spectra does not contain id column!")

  }

  if(!length(unique(x$id)) == length(x)) {

    stop("Spectra shall contain only unique entries!")

  }

  l <- list()

  for(method in methods) {

    adj_spec <- Spectra::compareSpectra(x,
                               FUN = get(method),
                               ...)

    colnames(adj_spec) <- x$id
    rownames(adj_spec) <- x$id

    ## add object to l
    new_index <- length(l) + 1
    l[[new_index]] <- adj_spec

    ## assign the name to the newly added entry
    names(l)[new_index] <- method

  }

  l

}
