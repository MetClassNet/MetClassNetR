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
  if(!"uniqueID" %in% spectraVariables(x)) {

    stop("Spectra does not contain uniqueID column!")

  }

  if(!length(x$uniqueID) == length(x)) {

    stop("Spectra shall contain only unique entries!")

  }

  l <- list()

  for(method in methods) {

    adj_spec <- compareSpectra(x,
                               FUN = get(method),
                               ...)

    colnames(adj_spec) <- x$uniqueID
    rownames(adj_spec) <- x$uniqueID

    ## add object to l
    new_index <- length(l) + 1
    l[[new_index]] <- adj_spec

    ## assign the name to the newly added entry
    names(l)[new_index] <- method

  }

  l

}
