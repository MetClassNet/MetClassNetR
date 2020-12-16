#' @name export2gml
#'
#' @aliases export2gml
#'
#' @title Create a network and export to gml from adjacency list
#'
#' @description
#' The function `export2gml` creates a network from adjacency lists.
#' The gml is saved in the current working directory having a
#' selected name (`name`). All or selected edge attributes may be
#' selected using the parameter `select`.
#'
#' @param x
#' `data.frame` adjacency list having `Var1` and `Var2` in the first
#' two columns and additional (edge) attributes in following columns.
#'
#' @param file
#' `character`, the name which the gml should have, e.g. "mz"
#'  produces the file 'mz.gml' in the current working directory.
#'
#' @param select
#' `character`, the default is `FALSE` and returns a gml having edge
#' attributes of all availabe columns of `x`. A `character` value or
#' string may be used to select one ore more columns by their columnnames
#' that should be saved as edge attributes.
#'
#' @details
#' ´`export2gml` uses adjacency matrixes to create first a network
#'  and then saving it as gml in the current working directory.
#'  The name of the gml file is selected by the User with the
#'  `file` parameter.
#'  The networks contains the feature ID's as nodes. The edges
#'  contains corresponding source  feature ID (sourceID) and target
#'  feature ID (targetID) and either all column attributes as additional
#'  edge attributes or selected ones if `select` was used.
#'
#' @return
#' `.gml` file saved in current working directory with selected name
#' using `file `.

#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#'  ######  example to be addded!
#'
#'
#'@export
exportNet2gml <- function (x, file, select = F, ...) {

  if (select != FALSE) {
    x <- x[,c("Var1", "Var2", select)]
  }

  net <- igraph::graph_from_data_frame(x, directed = TRUE, vertices = NULL)
  igraph::E(net)$sourceID <-  as.character(x$Var1)
  igraph::E(net)$targetID <-  as.character(x$Var2)
  fl <- paste(file, ".gml", sep="")
  igraph::write_graph(net, file = fl, format = c("gml"))
}


