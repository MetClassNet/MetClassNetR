
# ################################################################################
# # rename files
# RenameFiles <- function(netDir, patternToRemove = "") {
#   netDir <-
#     ifelse(
#       substr(netDir, nchar(netDir), nchar(netDir)) == "/",
#       substr(netDir, 1, (nchar(netDir)-1)),
#       netDir)
#
#   # list all files of the given format in the corresponding directory
#   files <-
#     list.files(
#       path = netDir,
#       full.names = TRUE,
#       pattern =  paste0(".*", patternToRemove, "(.*).gml")
#     )
#
#   # rename files
#   lapply(
#     files,
#     function(X) {
#       newName <- str_replace(X, paste0(".*", patternToRemove, "(.*.gml)"), "\\1")
#       file.rename(X, paste0(netDir, "/", newName))
#     }
#   )
#   return()
# }


# Function to generate 4 toy networks
# INPUT:
#   netDir - directory to store the toy networks
# OUTPUT: none, but it generates 4 files in netDir
#' @name makeToyNet
#'
#' @aliases makeToyNet
#'
#' @title Generate toy networks
#'
#' @description
#' The function `makeToyNet` generates 4 toy networks.
#'
#' @param netDir
#' `Path`, directory to store the toy networks
#'
#' @return
#' Nothing, but it creates four new files in `netDir`
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @examples
#' # See the NetworkComparison vignette
#'
#' @export
makeToyNet <- function(netDir) {
  # create networks
  n1 <-
    matrix(
      data = c("A", "B", "B", "C", "B", "E", "C", "D"),
      ncol = 2, byrow = TRUE)
  n2 <-
    matrix(
      data = c("A", "B", "B", "F", "A", "F", "C", "F", "C", "E", "C", "D"),
      ncol = 2, byrow = TRUE)
  n3 <-
    matrix(
      data = c("A", "F", "E", "F", "C", "F", "C", "E"), ncol = 2, byrow = TRUE)
  n4 <-
    matrix(
      data =
        c("A", "B", "A", "F", "B", "F", "B", "C", "A", "D", "C", "D",
          "F", "D"),
      ncol = 2, byrow = TRUE)

  # set column names
  colnames(n1) <-
    colnames(n2) <- colnames(n3) <- colnames(n4) <- c("node1", "node2")

  # create target directory
  dir.create(netDir)

  # save networks in files
  write.csv(
    n1, file = paste0(netDir, "Network1.csv"), row.names = FALSE,
    quote = FALSE)
  write.csv(
    n2, file = paste0(netDir, "Network2.csv"), row.names = FALSE,
    quote = FALSE)
  write.csv(
    n3, file = paste0(netDir, "Network3.csv"), row.names = FALSE,
    quote = FALSE)
  write.csv(
    n4, file = paste0(netDir, "Network4.csv"), row.names = FALSE,
    quote = FALSE)

  return()
}


# Function to read all the networks stored in a given directory
# INPUTS:
#   netDir   - directory containing all the networks.
#              NOTE. The networks must be stored in csv format
#   directed - boolean value indicating if the networks are directed or not,
#              the default value is FALSE (i.e., undirected networks)
#   pattern  - (optional) pattern that the file names must contain to be read.
#              If all the networks in netDir are to be read, pattern can be
#              omitted
#   format   - files' format. NOTE. If the format is "csv", then the files must
#              contain 2 columns (source - target) and each row must be an edge
# OUTPUT:
#   list of igraph objects, one per network to compare
#' @name readNet
#'
#' @aliases readNet
#'
#' @title Read all the networks stored in a given directory
#'
#' @description
#' The function `readNet` reads all the networks stored in a given directory
#'
#' @param netDir
#' `Path`, directory containing all the networks. NOTE. The networks must be
#' stored in csv format
#'
#' @param directed
#' `boolean`, value indicating if the networks are directed or not, the default
#' value is FALSE (i.e., undirected networks)
#'
#' @param  pattern
#' `string`, (optional) pattern that the file names must contain to be read.
#' NOTE. If all the networks in `netDir` are to be read, the parameter
#' `pattern` can be omitted
#'
#' @param  format
#' `string`, files' format (i.e., extension of the network files). NOTE. If the
#' format is "csv", then the files must contain 2 columns (source - target) and
#' each row must be an edge
#'
#' @import igraph
#'
#' @return
#' List of igraph objects, one per network to compare
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @examples
#' # See the NetworkComparison vignette
#'
#' @export
readNet <- function(netDir, directed = FALSE, pattern = "", format = "csv") {

  # if exists, remove final "/" from the network directory name
  netDir <-
    ifelse(
      substr(netDir, nchar(netDir), nchar(netDir)) == "/",
      substr(netDir, 1, (nchar(netDir)-1)),
      netDir)

  # list all files of the given format in the corresponding directory
  files <-
    list.files(
      path = netDir,
      pattern = paste0(".*", pattern, ".*.", format, "$"),
      full.names = TRUE
      )

  networks <- list()

  # loop through the files
  for (i in files) {
    # verify file format
    if (format == "csv") {
      # read file content
      n <- as.matrix(read.csv(i))

      # create igraph object
      net <- igraph::graph_from_edgelist(el = n, directed = directed)

    } else {
      net <- igraph::read_graph(i, format = format)

      if (directed == FALSE) {
        net <- igraph::as.undirected(net)
      }
    }

    # add network to named list
    networks[[str_replace(i, paste0(".+/(.*).", format), "\\1")]] <- net
  }

  return(networks)
}


# Function to calculate, create a table and plot the following statistics:
#    - Density (no self loops are considered)
#    - Diameter (if there is more than one connected component, the largest
#        diameter will be returned, independently of the size of the connected
#        component)
#    - Average degree
#    - Average path length
#    - Clustering coefficient (the networks are considered as undirected)
# It also plots the degree distribution, and upset plots of the overlap of
# nodes and edges
# INPUT:
#   net - list of igraph objects, as returned by readNet()
# OUTPUT:
#   data frame with the statistics
#' @name calculNetStats
#'
#' @aliases calculNetStats
#'
#' @title Calculate, create a table and plot statistics
#'
#' @description
#' The function `calculNetStats` serves to calculate, create a table and plot
#' the following statistics:
#'    - Density (no self loops are considered)
#'    - Diameter (if there is more than one connected component, the largest
#'      diameter will be returned, independently of the size of the connected
#'      component)
#'    - Average degree
#'    - Average path length
#'    - Clustering coefficient (the networks are considered as undirected)
#' It also plots the degree distribution, and upset plots of the overlap of
#' nodes and edges reads all the networks stored in a given directory
#'
#' @param net
#' `list`, igraph objects, as returned by readNet()
#'
#' @import igraph
#'
#' @return
#' Data frame with the statistics
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @examples
#' # See the NetworkComparison vignette
#'
#' @export
calculNetStats <- function(net) {

  # initialize variables
  netStats <-
    data.frame(name = character(), vcount = double(), ecount = double(),
      dens = double(), diameter = double(), avDeg = double(), noCC = double(),
      avPathLen = double(), clustCo = double())

  netClo <-
    data.frame(network = character(), node = character(), closeness = double())

  netBet <-
    data.frame(network = character(), node = character(),
      betweenness = double())

  netCC <- list()

  BCCStats <-
    data.frame(name = character(), vcount = double(), ecount = double(),
      dens = double(), diameter = double(), avDeg = double(), noCC = double(),
      avPathLen = double(), clustCo = double())

  BCCClo <-
    data.frame(network = character(), node = character(), closeness = double())

  BCCBet <-
    data.frame(network = character(), node = character(),
      betweenness = double())

  noNet <- length(net) # number of networks

  netDeg <- getNodeDeg(net) # get nodes' degrees of all networks

  for (i in seq_len(noNet)) { # loop through networks to calculate stats
    # -- stats of the whole network
    netName <- names(net)[i] # get network's name
    n <- net[[i]] # get network
    cc <- components(n) # get connected components (CCs)
    ccSize <- makeTableCCSize(cc) # table of CCs sizes
    netCC[[netName]] <- cc # save CCs
    netStats <- rbind(netStats, basicNetStats(n, netName)) # basic stats
    netClo <- rbind(netClo, getClosenessCC(cc, n, netName)) # closeness
    netBet <- rbind(netBet, getBetweenness(n, netName)) # betweenness

    # -- stats of the biggest connected component (BCC)
    # get BCC
    BCC <-
      igraph::induced_subgraph(n, V(n)[cc$membership == which.max(cc$csize)])
    BCCStats <-
      rbind(BCCStats, basicNetStats(BCC, paste0(netName, "_BCC"))) # stats
    BCCClo <-
      rbind(BCCClo, getClosenessCC(components(BCC), n, netName)) # closeness
    BCCBet <- rbind(BCCBet, getBetweenness(BCC, netName)) #  betweenness
  }

  allStats <-
    list(netStats = netStats, networksDegrees = netDeg, networkCC = netCC,
      ccSize = ccSize, networkCloseness = netClo, networkBetweenness = netBet,
      BCC = BCC, BCCStats = BCCStats, BCCClo = BCCClo, BCCBet = BCCBet)

  return(allStats)
}


# Function to calculate basic network stats (density, diameter, average degree,
# average path length, and clustering coefficient)
# INPUTS:
#   network - igraph object
#   netName - name of the network
# OUTPUT: data frame with all the stats
#' @name basicNetStats
#'
#' @title Calculate basic network stats
#'
#' @import igraph
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
basicNetStats <- function(network, netName) {
  # fill new data frame with stats
  netStats <-
    data.frame(
      name = netName,
      vcount = igraph::vcount(network),
      ecount = igraph::ecount(network),
      dens = igraph::edge_density(network),
      diameter = igraph::diameter(network, unconnected = TRUE),
      avDeg = mean(degree(network)),
      noCC = igraph::components(network)$no,
      avPathLen =
        igraph::mean_distance(network, directed = is_directed(network)),
      clustCo = igraph::transitivity(network, type = "global")
      )

  return(netStats)
}


# Function to get the nodes' degrees
# INPUT:
#   net - list of networks as igraph objects
# OUTPUT: data frame with the nodes' degrees of all the networks in the list
#' @name getNodeDeg
#'
#' @title Get the nodes' degrees
#'
#' @import igraph
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
getNodeDeg <- function(net) {

  netDeg <-
    data.frame(network = character(), node = character(), degree = double())

  # loop through all the networks
  for (i in seq_len(length(net))) {
    # get nodes' degrees
    deg <-
      data.frame(
        network = names(net)[i],
        node = names(igraph::V(net[[i]])),
        degree = igraph::degree(net[[i]])
      )
    netDeg <- rbind(netDeg, deg)
  }

  return(netDeg)
}


# Function to make a table of the sizes of the connected components
# INPUT:
#    cc - connected components
# OUTPUT: table of sizes
#' @name makeTableCCSize
#'
#' @title Make a table of the sizes of the connected components
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
makeTableCCSize <- function(cc) {

  # make a table to check the size of the connected components
  ccSize <- as.data.frame(table(cc$csize))

  # rename column
  colnames(ccSize)[1] <- "ComponentSize"

  # leave only components with at least 2 nodes
  ccSize <- ccSize[as.numeric(as.character(ccSize$ComponentSize)) > 1,]

  return(ccSize)
}


# Function to calculate the closeness of each connected component of a network
# INPUTS:
#   cc      - connected components
#   network - network as igraph object
#   netName - network name
# OUTPUT: closseness of all connected components
#' @name getClosenessCC
#'
#' @title Calculate the closeness of each connected component of a network
#'
#' @import igraph
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
getClosenessCC <- function(cc, network, netName) {

  # calculate the closeness in each component
  clo <-
    unlist(
      lapply(
        seq_len(cc$no),
        function(X) {

          # get nodes' IDs of the nodes in current connected component
          nodesIDs <-
            igraph::V(network)[which(names(igraph::V(network)) %in%
              names(cc$membership[cc$membership == X]))]

          # calculate closeness of current Connected Component (CC)
          # NOTE. If normalized, the most central node (i.e., node with the
          # shortest path to the rest of the nodes) has the highest closeness,
          # or the lowest value otherwise
          # NOTE 2. Each CC is considered as an independent graph
          igraph::closeness(
            igraph::induced_subgraph(network, vids = nodesIDs),
            normalized = TRUE
            )
        }
      )
    )

  # add closeness data
  netClo <- data.frame(network = netName, node = names(clo), closeness = clo)

  return(netClo)
}


# Function to calculate the betweenness of the nodes in a network
# INPUTS:
#   cc      - connected components
#   network - network as igraph object
#   netName - network name
# OUTPUT: data frame containing the betweenness of the nodes in the network
#' @name getBetweenness
#'
#' @title Calculate the betweenness of the nodes in a network
#'
#' @import igraph
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
getBetweenness <- function(network, netName) {

  # calculate the betweenness
  bet <-
    data.frame(
      network = netName,
      node = names(igraph::V(network)),
      betweenness = igraph::betweenness(network)
      )

  return(bet)
}


# Function to make and print several plots
# INPUT:
#   stats - list of results, as returned by the calculNetStats function
# OUTPUT:
#   none, but it prints several plots
#' @name printStatsPlots
#'
#' @aliases printStatsPlots
#'
#' @title Calculate, create a table and plot statistics
#'
#' @description
#' The function `printStatsPlots` makes and prints several plots
#'
#' @param stats
#' `list`, results, as returned by the `calculNetStats` function
#'
#' @return
#' Nothing, but it prints several plots
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @examples
#' # See the NetworkComparison vignette
#'
#' @export
printStatsPlots <- function(stats) {
  # get data
  netStats <- stats[['netStats']]
  netDeg <- stats[['networksDegrees']]
  netClo <- stats[['networkCloseness']]
  netBet <- stats[['networkBetweenness']]
  ccSize <- stats[['ccSize']]
  BCC <- stats[['BCC']]
  BCCStats <- stats[['BCCStats']]
  BCCClo <- stats[['BCCClo']]
  BCCBet <- stats[['BCCBet']]
  isolNodes <- as.data.frame(table(netDeg[netDeg$degree == 0, c(1, 3)]))
  colnames(isolNodes)[1] <- "name"

  # make scatter plots with some of the data
  makeScatterPlot(isolNodes, stat2Plot = "Freq",
    title = "Number of isolated nodes", printLab = TRUE)

  s2p <- c("noCC", "dens", "diameter", "avPathLen", "clustCo", "avDeg")
  t <- c("Number of connected components", "Density", "Diameter",
    "Average path length", "Clustering coefficient", "Average degree")
  p <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)

  for (i in seq_len(length(s2p))) {
    makeScatterPlot(
      netStats, stat2Plot = s2p[i], title = t[i], printLab = p[i])
  }

  makeBoxPlot(netDeg, stat2Plot = "degree", log = TRUE)
  makeHist(
    netDeg, stat2Plot = "degree", binWidth = 1, minDeg = 0,
    maxDeg = max(netDeg$degree))
  makeBoxPlot(
    netClo, stat2Plot = "closeness", bestValue = "max", printLab = TRUE)
  makeBoxPlot(
    netBet, stat2Plot = "betweenness", bestValue = "max", printLab = TRUE)

  return()
}


# Function to make and print a scatter plot
# INPUTS:
#   netStats  - data frame, as returned as part of the result of
#               calculNetStats
#   stat2Plot - string defining the stat that wishes to be plotted. It can be
#               "density", "diameter", "avPathLength" or "clustCo"
# OUTPUT:
#   none, but it prints the corresponding scatter plot
#' @name makeScatterPlot
#'
#' @title Make and print a scatter plot
#'
#' @import ggplot2
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
makeScatterPlot <- function(netStats, stat2Plot, title ="", printLab = FALSE) {
  printNothing(1)
  p <-
    ggplot2::ggplot(
      netStats,
      ggplot2::aes(
        x = name,
        y = eval(as.name(stat2Plot)),
        color = rainbow(nrow(netStats))
        )
      ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::labs(
      title =
        ifelse(
          title != "",
          paste0(title, " comparison\n"),
          paste0(stringr::str_to_title(stat2Plot), " comparison\n")
          ),
      x = "Network",
      y = ifelse(title != "", title, stringr::str_to_title(stat2Plot))
      ) +
    ggplot2::theme(
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  # check if labels are to be printed
  if (printLab == TRUE) {
    # check if the values to plot contain floating point
    if (typeof(netStats[[stat2Plot]]) == "double") {
      # keep only 6 decimals for the labels
      labels <-
        sapply(
          netStats[[stat2Plot]],
          function(X) {
            format(round(X, 6), nsmall = 5)
          }
        )
    } else {
      labels <- netStats[[stat2Plot]]
    }
    p <- p + ggplot2::geom_text(label = labels, hjust = -0.3)  # add labels
  }
  print(p) # print plot

  return()
}



# Function to make and print a boxplot
# INPUTS:
#   data      - data frame, as returned as part of the result of
#               calculNetStats
#   stat2Plot - string defining the stat that wishes to be plotted. It can be
#               "degree" or "closeness"
#   bestValue - string to define whether the highest or lowest values are the
#               best ones. This is only useful if printLab == TRUE. Possible
#               values are "max" and "min"
#   printLab -  flag to choose whether to show the label of the highest/lowest
#               value
# OUTPUT:
#   none, but it prints the corresponding scatter plot
#' @name makeBoxPlot
#'
#' @title Make and print a boxplot
#'
#' @import ggplot2
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
makeBoxPlot <- function(data, stat2Plot, bestValue = "max", printLab = FALSE,
  log = FALSE) {

  printNothing(2)

  # remove NaN's
  data <- data[!is.nan(data[[stat2Plot]]), ]

  # check if logarithmic transformation is to be made
  if (log == TRUE) {
    data[[stat2Plot]] <- log(data[[stat2Plot]] + 1)
  }

  p <-
    ggplot2::ggplot(
      data,
      ggplot2::aes(x = network, y = eval(as.name(stat2Plot)), fill = network)
      ) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      title = paste0(stringr::str_to_title(stat2Plot), " comparison\n"),
      x = "Network",
      y = stringr::str_to_title(stat2Plot)
      ) +
    ggplot2::theme(
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  # verify if labels should be added to the plot
  if(printLab == TRUE) {
    plotLabels <- getPlotLabels(data, stat2Plot, bestValue)
    p <- p + ggplot2::geom_text(data = plotLabels, aes(label = node),
      position = "identity", vjust = ifelse(bestValue == "max", -0.3, 0.3),
      size = 2)
  }
  print(p)

  return()
}



# Function to get the labels to plot in the box plot
# INPUTS:
#   data      - data frame obtained via the calculNetStats function
#   stat2Plot - string defining the stat that wishes to be plotted. It can be
#               "degree" or "closeness"
#   bestValue - string to define whether the highest or lowest values are the
#               best ones. This is only useful if printLab == TRUE. Possible
#               values are "max" and "min"
#' @name getPlotLabels
#'
#' @title Get the labels to plot in the box plot
#'
#' @importFrom plyr ddply
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
getPlotLabels <- function(data, stat2Plot, bestValue) {
  plotLabels <-
    data[as.logical(ave(data[, eval(stat2Plot)], data$network,
    FUN = function(x) x == eval(str2expression(paste0(bestValue, "(x)"))))), ]

  # check if there are more than 2 nodes' names to be printed for any network
  if (any(table(plotLabels$network) > 2)) {
    newLab <- data.frame() # initialize empty data frame

    for (net in unique(plotLabels$network)) { # loop through the network names
      rows <- which(plotLabels$network == net) # get rows from current network

      # print some of the labels for the current network
      print(
        paste0("Network ", net, " has ", length(rows), " nodes with ",
        bestValue, " ", stat2Plot, ". Here are some of them: ")
        )
      print(plotLabels$node[head(rows)])
      printNothing(1)

      if (length(rows) > 2) { # if there are more than 2 rows
        selRows <- sample(rows, 2)  # pick two rows at random
        newLab <- rbind(newLab, plotLabels[selRows,]) # add labels
      } else {
        newLab <- rbind(newLab, plotLabels[rows,])
      }
    }
    plotLabels <- newLab
  }

  # put all the nodes' names to be printed per network together
  plotLabels <-
    plyr::ddply(
      plotLabels,
      .(network),
      function(x) {
        return(
          c(node = paste(x$node, collapse = ", "), x = unique(x[[stat2Plot]]))
          )
        }
      )
  colnames(plotLabels)[ncol(plotLabels)] <- stat2Plot
  plotLabels[[stat2Plot]] <- as.numeric(plotLabels[[stat2Plot]] )

  return(plotLabels)
}


# Function to print empty lines
# INPUT: number of empty lines to print
# OUTPUT: none but prints the lines
#' @name printNothing
#'
#' @title Print empty lines
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
printNothing <- function(n) {
  for (i in seq_len(n)) {
    print("")
  }

  return()
}

# Function to calculate and print a set of histograms
# INPUTS:
#   data      - data frame, as returned as part of the result of
#               calculNetStats
#   stat2Plot - string defining the stat that wishes to be plotted. It can be
#               "degree"
#   binWidth  - width of the bins for each histogram. Default value = 1
#   minDeg    - minimum degree to consider for the plot. Default value = 0
#   maxDeg    - maximum degree to consider for the plot. Default value = 10
# OUTPUT:
#   none, but it prints the corresponding histograms
#' @name makeHist
#'
#' @title Calculate and print a set of histograms
#'
#' @import ggplot2
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
makeHist <- function(data, stat2Plot, binWidth = 1, minDeg = 0, maxDeg = 10) {

  # filter data
  data <-
    data %>%
    filter(
      eval(as.name(stat2Plot)) >= minDeg & eval(as.name(stat2Plot)) <= maxDeg
      )

  printNothing(1)
  print(
    ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = eval(as.name(stat2Plot)), color = network, fill = network)
    ) +
    ggplot2::geom_histogram(
      alpha = 0.5, binwidth = binWidth, position = "identity")
  )

  return()
}


# Function to calculate and print plots of the overlap of nodes and edges
# between a set of networks
# INPUT:
#   networks - list of igraph objects to analyze
# OUTPUT:
#   none, but it prints the corresponding plots
#' @name calculateOverlap
#'
#' @aliases calculateOverlap
#'
#' @title Calculate, create a table and plot statistics
#'
#' @description
#' The function `calculateOverlap` serves to calculate and print plots of the
#' overlap of nodes and edges
#'
#' @param networks
#' `list`, igraph objects to analyze
#'
#' @import igraph
#'
#' @return
#' Nothing, but it prints the corresponding plots
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @examples
#' # See the NetworkComparison vignette
#'
#' @export
calculateOverlap <- function(networks) {
  # get the number of networks
  noNet <- length(networks)

  # initialize empty list to save the nodes/edges of each network
  allNodes <- list()
  allEdges <- list()

  # loop through the networks to get the nodes' and edges' list
  for(i in seq_len(length(networks))) {
    # save nodes' names
    allNodes[[names(networks)[i]]] <- names(igraph::V(networks[[i]]))

    # get edges' list
    edges <- igraph::as_edgelist(networks[[i]])

    # sort edges alphabetically
    edges <- t(apply(edges, 1, sort))

    # paste the edges
    edges <- apply(edges, 1, paste, collapse = ".")

    # save edges
    allEdges[[names(networks)[i]]] <- edges
  }

  printNothing(1)

  # make upset plot for the overlap of nodes
  makeUpsetPlot(
    allNodes, title = "Overlap of nodes",
    xLabel = "Nodes per network", yLabel = "Nodes' intersection"
    )

  printNothing(1)

  # make upset plot for the overlap of edges
  makeUpsetPlot(
    allEdges, title = "Overlap of edges",
    xLabel = "Edges per network", yLabel = "Edges' intersection"
    )

  return()
}



# Function to get the overlapping nodes from a given set of input networks
# INPUTS:
#   networks      - list of networks (igraph objects)
#   networksIndex - list of networks' index to consider for the overlap, e.g.,
#                   c(1, 3) to take the first and third networks
# OUTPUT:
#   list of overlapping nodes
#' @name getOverlappingNodes
#'
#' @aliases getOverlappingNodes
#'
#' @title Get the overlapping nodes
#'
#' @description
#' The function `getOverlappingNodes` gets the overlapping nodes from a given
#' set of input networks
#'
#' @param networks
#' `list`, list of networks (igraph objects)
#'
#' @param networksIndex
#' `list`, list of networks' index to consider for the overlap, e.g.,
#' c(1, 3) to take the first and third networks
#'
#' @import igraph
#'
#' @return
#' List of overlapping nodes
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @examples
#' # See the NetworkComparison vignette
#'
#' @export
getOverlappingNodes <- function(networks, networksIndex) {
  if (exists("networksIndex")) {
    # filter network list to keep only the networks to process
    networks <- networks[networksIndex]
  }

  # get nodes' names
  nodes <- lapply(networks, function(X) {names(igraph::V(X))})

  # obtain the intersection
  overlappingNodes <- Reduce(intersect, nodes)

  return(overlappingNodes)
}


# Function to make and print an upset plot
# INPUTS:
#   dataToPlot - named list containing the data to plot
#   xLabel     - string defining the label for the X axis, per default it is
#                "Set Size"
#   yLabel -     string defining the label for the Y axis, per default it is
#                "Intersection Size"
# OUTPUT:
#   none, but it prints the corresponding upset plot
#' @name makeUpsetPlot
#'
#' @title Make and print an upset plot
#'
#' @importFrom UpSetR upset fromList
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
#' @export
makeUpsetPlot <- function(dataToPlot, title = "Overlap", xLabel = "Set Size",
  yLabel = "Intersection Size") {

  printNothing(1)

  print(
    UpSetR::upset(
      UpSetR::fromList(dataToPlot), order.by = "freq", point.size = 2,
      line.size = 1, mainbar.y.label = yLabel, sets.x.label = xLabel,
      #empty.intersections = "on",
      # y-axis label, y-axis ruler numbers, x-axis label, x-axis ruler numbers,
      # x-axis legend content, overlap size legend
      text.scale = c(1.3, 0.9, 1, 1, 1.1, 0.9)
    )
  )

  grid::grid.text(title, x = 0.65, y = 0.95, gp = grid::gpar(fontsize = 16))

  return()
}
