# remove compartment from GSMN nodes' labels
g <- multiLayer$layers[["gsmn"]]
lab <- vertex_attr(g, "label")
newLab <- sapply(lab, function(X) {substr(X, 1, nchar(X) - 2)})
g <- set_vertex_attr(g, "label", index = V(g), newLab)
g <- set.vertex.attribute(g, "name", index = V(g), newLab)


# set new node names
multiLayer$layers[["gsmn"]] <- g

allEdges <-
  data.frame(
    node1 = character(),
    node2 = character(),
    source = character()
    )

# get edge list of all the networks and paste it in a single data frame
for (i in seq_len(length(multiLayer$layers))) {
  # get edge list
  edges <- as_edgelist(multiLayer$layers[[i]])

  # add to general list
  allEdges <-
    rbind(
      allEdges,
      data.frame(
        node1 = edges[, 1],
        node2 = edges[, 2],
        source =
          paste(multiLayer$type[i], names(multiLayer$layers)[i], sep = "_")
      )
    )
}

# add inter-layer edges
allEdges <-
  rbind(
    allEdges,
    data.frame(
      node1 = multiLayer$interLayerEdges[, 1],
      node2 = multiLayer$interLayerEdges[, 2],
      source = "InterLayer"
    )
  )

write.csv(allEdges, "~/Documents/MetClassNet/R/MetClassNetR_devel_Liesa_Elva/MultiLayer4Cytoscape.csv", row.names = FALSE, quote = FALSE)

expNodes <- unique(unlist(lapply(multiLayer$layers[multiLayer$type == "Exp"], function(X) {names(V(X)) })))
gsmnNodes <- unique(unlist(lapply(multiLayer$layers[multiLayer$type == "GSMN"], function(X) {names(V(X)) })))
allNodes <- data.frame(node = c(expNodes, gsmnNodes), nodeType = c(rep("Exp", length(expNodes)), rep("GSMN", length(gsmnNodes))))


write.csv(allNodes, "~/Documents/MetClassNet/R/MetClassNetR_devel_Liesa_Elva/MultiLayer4Cytoscape_NodeList.csv", row.names = FALSE, quote = FALSE)


