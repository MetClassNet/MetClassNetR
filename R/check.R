#'
#'
#' @export
checkQFeatures <- function(x) {

  if(length(rowDataNames(mtbls1586_qf_neg)) > 1) {

    stop("Currently only QFeatures with a single assay are supported")

  }

  if(!all(c("rtime", "mz", "id") %in% rowDataNames(x)[[1]])) {

    stop("Missing minimal column header: 'rtime', 'mz', 'id'")

  }

  TRUE

}

#'
#'
#' @export
checkSpectra <- function(x) {

  # check if id is present
  if(!c("id") %in% spectraVariables(x)) {

    stop("Missing minimal metadata: 'id'")

  }

  # check if only a single spectra per id exist
  if(!length(unique(x$id)) == length(x)) {

    stop("Only single MS2 spectra per feature allowed, perform consolidation of MS2 first!")

  }

  TRUE

}

#'
#'
#' @export
checkIdNamespace <- function(QFeatures, Spectra) {

  # get IDs from QFeatures
  features_ids <- sort(unique(rownames(QFeatures)[[1]]))

  # get IDs from Spectra
  spectra_ids <- sort(unique(Spectra$id))

  if(!length(features_ids) == length(spectra_ids) && !all(features_ids == spectra_ids)) {

    stop("Ids between QFeatures and Spectra are not matching")

  }

  TRUE

}
