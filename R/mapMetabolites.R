####
##
#  To Do
#
#   Additional parameter descriptions in the documentation of the main function oder of the single functions which are not exported?
#   An der Stelle die Beschreibungen der zusätzlichen Parameter hinzufügen die bei ... eingesetzt werden können. Entweder hier als
#   zusätzliche Erklärung. Oder prüfen ob die Doku zu den nicht exportiereten Funktionen trotzdem eingesehen werden kann.
#
#  examples section either fill with code or remove. In addition to code the sentecne detailed example in vignette might fit.
#
#  package check errors
#  package build errors
#
#
#  TEST TEST TEST
#
#
###


#' @name mapMetToGSMN
#'
#' @aliases mapMetToGSMN
#'
#' @title Map the experimental nodes to GSMN nodes
#'
#' @description
#' Function to map the identified experimental nodes to the
#' corresponding GSMN nodes.
#'
#' @param inputData
#' `list`, list returned by the `loadInputData` function
#'
#' @param method one of MetExplore, molecularFormula, accurateMass, class_chebi or class_chemont
#'
#' @param ... additional parameter depending on the chosen function
#'
#' @return
#' Data frame with the mappings and ontology-based distances.
#'
#' @author Sarah Scharfenberg
#'
#' @export
mapMetToGSMN <- function(inputData, ...) {

  if(method=="metexplore"){

    ## check parameter
    ## is it possible to check the existence of an paramter in ...
    #if(resFile)

    ## call function
    mapMetToGSMN_metexplore(inputData, ...)
  }

  if(method=="molecularFormula"){

    ## check parameter
    ## call function
    mapMetToGSMN_molecularFormula()
  }

  if(method=="accurateMass"){

    ## check parameter
    ## call function

    ## Michael add code here
    mapMetToGSMN_accurateMass()
  }

  if(method=="class"){

    ## check parameter
      ## obo
      ## take metabolites form input or input full list?

    ## call function
      ## one per chebi and chemont??
      ## update this in documentation above
    mapMetToGSMN_class()
  }
}




#' @name mapMetToGSMN_metexplore
#'
#' @aliases mapMetToGSMN
#'
#' @title Map the experimental nodes to GSMN nodes
#'
#' @description
#' Function to map the identified experimental nodes to the
#' corresponding GSMN nodes, using the ChEBI ontology. This function calls the
#' Metabolomics2Networks Python package
#'
#' @param inputData
#' `list`, list returned by the `loadInputData` function
#'
#' @param resFile
#' `character`, file name for the resulting mappings file. The default value is
#' "Res_Met2Net_MappedMet.txt"
#'
#' @return
#' Data frame with the mappings and ontology-based distances.
#'
#' @author Elva Novoa, \email{elva-maria.novoa-del-toro@@inrae.fr}
#'
mapMetToGSMN_metexplore <- function(inputData, resFile = "Res_Met2Net_MappedMet.txt") {

  # generate command line to execute metabolites2Network
  com <-
    paste0(
      "python3 ",
      inputData$met2NetDir,
      "metabolomics2network.py tsv ",
      inputData$idenMetF,
      " ",
      inputData$compF,
      " ",
      inputData$resPath,
      resFile,
      " ",
      inputData$configF,
      " 2"
    )

  # execute code
  system(com)

  return()
}
