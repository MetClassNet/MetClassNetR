#' @name exportNet2gml
#'
#' @aliases exportNet2gml
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
#' `data.frame` adjacency list having `Row` and `Col` in the first
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
#' ´`exportNet2gml` uses adjacency matrices to create first a network
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
#'  ######  example to be added!
#'
#'
#'@export
exportNet2gml <- function (x, file, select = F, ...) {

  if (select != FALSE) {
    x <- x[,c("Row", "Col", select)]
  }

  net <- igraph::graph_from_data_frame(x, directed = TRUE, vertices = NULL)
  igraph::E(net)$sourceName <-  as.character(x$Row)
  igraph::E(net)$targetName <-  as.character(x$Col)
  fl <- paste(file, ".gml", sep="")
  igraph::write_graph(net, file = fl, format = c("gml"))
}

#' @name exportAttributes2gml
#'
#' @aliases exportAttributes2gml
#'
#' @title Create a network and export to gml from adjacency list
#'
#' @description
#' The function `exportAttributes2gml` creates a network from adjacency lists.
#' The gml is saved in the current working directory having a
#' selected name (`name`). All or selected edge attributes may be
#' selected using the parameter `select`.
#'
#' @param x
#' `data.frame` adjacency list having `Row` and `Col` in the first
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
#' ´`exportNet2gml` uses adjacency matrices to create first a network
#'  and then saving it as gml in the current working directory.
#'  The name of the gml file is selected by the User with the
#'  `file` parameter.
#'  The networks contains the feature ID's as nodes. The edges
#'  contains corresponding source  feature ID (sourceID) and target
#'  feature ID (targetID) and either all column attributes as additional
#'  edge attributes or selected ones if `select` was used.
#'  `names` contains RT, m/z and manual Annotations that will be stored as node attributes
#'
#' @return
#' `.gml` file saved in current working directory with selected name
#' using `file `.

#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#'  ######  example to be added!
#'
#'
#'@export
exportAttributes2gml <- function (x, file, select = F, anno, ...) {

  if (select != FALSE) {
    x <- x[,c("Row", "Col", select)]
  }

  net <- igraph::graph_from_data_frame(x, directed = TRUE, vertices = row.names(anno))
  igraph::E(net)$sourceName <-  as.character(x$Row)
  igraph::E(net)$targetName <-  as.character(x$Col)
  net <-  set_vertex_attr(graph = net, name = "mz", value = as.character(anno$mz))
  net <-  set_vertex_attr(graph = net, name = "rt", value = as.character(anno$RT))
  net <-  set_vertex_attr(graph = net, name = "Annotation", value = as.character(anno$Annotation))
  # net <- igraph::set_vertex_attr(net, "mz") <- as.character(x$Var1_mz)
  fl <- paste(file, ".gml", sep="")
  igraph::write_graph(net, file = fl, format = c("gml"))
}


