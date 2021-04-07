threshold <- function(statistical, type, args,
                      values = c("all", "min", "max"), ...) {
  ## checks if p values are included in list and deletes them
  if ( TRUE == any(endsWith(names(statistical), "_p"))) {
    l <- statistical[!endsWith(names(statistical), "_p")]
  }
  else {l <- statistical}
  ## args, either N for tops
  ## args, either N for tops
  ## or a list of threshold
  if (any(duplicated(names(args)))) {
    stop("names(args) contain duplicated entries")
  }
  ## check match.arg for values
  values <- match.arg(values)
  if (type == "threshold") {
    ## iterate through the list and remove the links below or above the
    ## threshold and write to list
    l <- lapply(seq_along(l), function(x) {
      ## find corresponding model in l
      name_x <- names(l)[x]
      ## get corresponding threshold in args
      threshold_x <- args[[names(l)[x]]]
      ## get corresponding adjacency matrix in l
      ## corresponds to MetNet function
      l_x <- l[[name_x]]
      ## for pearson/spearman correlation models (incl. partial and
      ## semi-partial), lasso, randomForest, clr, aracne and bayes higher
      ## values corresond to higher confidence
      ## only assign 1 to values that are above the threshold
      ifelse(abs(l_x) > threshold_x, l_x, NA)
    })
    ## allow for compatibility of arguments
    ## calculate consenses from the binary matrices
    cons <- MetNet:::threeDotsCall(sna::consensus, dat = l, ...)
    ## threshold consensus that it is a binary matrix
    cons <- ifelse(cons >= args$threshold, cons, NA)
    rownames(cons) <- colnames(cons) <- colnames(l[[1]])
  }
  names(l) <- names(statistical[!endsWith(names(statistical), "_p")])
  l[["Consensus"]] <- cons
  #class(l[[3]]) <- "numeric"
  return(l)
}
