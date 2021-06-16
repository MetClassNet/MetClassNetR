#' @name spec_adjacency_list
#'
#' @aliases spec_df
#'
#' @title Create a `data.frame` from spectral similarity matrix
#'
#' @description
#' The function `spec_df` creates a `data.frame` from
#' spectral similarity matrix
#'
#' @param x
#' `data.frame` adjacency matrix
#'
#'
#' @return
#' `data.frame` adjacency list
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#'  ######  example to be added!
#'
#'
#' @export
spec_df <- function(x){

  x[[1]][upper.tri(x[[1]])] <- ''

  simil <- reshape2::melt(x[[1]]) %>%
    filter(Var1 != Var2) %>%
    filter(value != '') %>%
    filter(value != "NaN" )


  colnames(simil) <- c("Row", "Col", "value")

  return(simil)
}
