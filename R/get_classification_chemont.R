get_classification_chemont <- function(network_tibble, delay = 0){

  net_classes <- sapply(network_tibble$Identifiers.inchikey,function(ik){
    Sys.sleep(delay)
    classific <- eval(classyfireR::get_classification(ik))
    if(!is.null(classific)){return(classific@classification)}
    return(NULL)
  })

  ## if all are classified correctly the output of sapply is not a list of tibbles
  ## but a hughe matrix with level, classification and chemont as lists of characters
  network_tibble <- network_tibble %>%
    add_column(classes = net_classes)

  # add 2nd column: find NULL entries and mark them as unclassified
  tmp <- sapply(network_tibble$classes,
                function(x){
                  return(!is.null(unlist(x)))
                  })
  network_tibble <- network_tibble %>%
    mutate(classified = tmp)
  rm(tmp)

  ## add 3rd column: report the most detailed class
  tmp <- unlist(sapply(network_tibble$classes,
                function(x){
                  if(!is.null(unlist(x))){
                    return(as.character(tail(x,n=1)$CHEMONT))
                  }
                  else{
                    return("")
                  }
                }))
  network_tibble <- network_tibble %>%
    mutate(class = tmp)
  rm(tmp)

  return(network_tibble)
}
