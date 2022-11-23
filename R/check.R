#'
#'
#' @export
checkQFeatures <- function(x) {

  if(length(rowDataNames(x)) > 1) {
    stop("Currently only QFeatures with a single assay are supported")
  }

  if(!all(c("rtime", "mz", "id") %in% rowDataNames(x)[[1]])) {
    stop("Missing minimal column header: 'rtime', 'mz', 'id'")
  }

  return(TRUE)

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

  return(TRUE)
}

#'
#'
#' @export
checkIdNamespace <- function(QFeatures, Spectra) {

  # get IDs from QFeatures
  features_ids <- unique(rownames(QFeatures)[[1]])

  # get IDs from Spectra
  spectra_ids <- unique(Spectra$id)

  if(!all(spectra_ids %in% features_ids)) {
    stop("Mismatch of features between the peak list file and the spectra one")
  }

  return(TRUE)
}
