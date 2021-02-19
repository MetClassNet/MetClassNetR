#' @name spec_adjacency_list
#'
#' @aliases spec_adjacency_list
#'
#' @title Create an adjacency list from spectral similarity matrix
#'
#' @description
#' The function `spec_adjacency_list` creates an adjacency lists from
#' spectral similarity matrix and is comparable to MetNet function adjacency_list()
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
spec_adjacency_list <- function(x){

  x[[1]][upper.tri(x[[1]])] <- ''

  simil <- reshape2::melt(x[[1]]) %>%
    filter(Var1 != Var2) %>%
    filter(value != '') %>%
    filter(value != "NaN" )


  colnames(simil) <- c("Var1", "Var2", "value")

  return(simil)
}
