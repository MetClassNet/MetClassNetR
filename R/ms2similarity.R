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
