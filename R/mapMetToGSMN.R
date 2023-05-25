###################################################
#      General call of the function               #
###################################################

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
#' @param method one of Metabolomics2Network, molecularFormula, accurateMass, class_chebi or class_chemont
#'
#' @param ... additional parameter depending on the chosen function
#'
#' @return
#' Data frame with the mappings and ontology-based distances.
#'
#' @author Elva Novoa, Sarah Scharfenberg
#'
#' @export
mapMetToGSMN <- function(inputData, method="metabolomics2network", resFile = "Res_Met2Net_MappedMet.txt", ...) {

  if(method=="metabolomics2network"){

    ## call function
    .mapMetToGSMN_metabolomics2network(inputData, resFile)
  }

  if(method=="id_inchikey"){

    ## check parameter


    ## call function
    .mapMetToGSMN_inchikey(inputData)
  }


}


###################################################
#      Single functions - not exported            #
###################################################


# The Python tool `Metabolomics2Network` (see
# https://forgemia.inra.fr/metexplore/metabolomics2network) uses the ChEBI
# ontology to localize a ChEBI id from the list of identified compounds, and it
# gives as result the ChEBI id of the nearest metabolite that is present in the
# GSMN, as well as the distance between them.
# In order to run `Metabolomics2Network`, we need to call the
# `mapMetToGSMN` function and pass the `inputData` object as input
# parameter. We can also choose the name of the resulting mappings file
# (`resFile`), which per default is "Res_Met2Net_MappedMet.txt".
# Such restulting file contains 12 columns, but we will focus on the first 5:
#   1. `metabolite name`. Name of our metabolites, as defined in the list of
# identified metabolites (`id` column from the file
#                         `IdentifiedMet_MTBLS1586_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv`)
# 2. `mapped on id`. Id of the node from the GSMN to which the metabolite was
# mapped (`ID` column from the file `WormJamMet.tsv`).
# 3. `mapping type`. Either `exact multimapping` or `chebi class mapping`,
# depending on which type we selected when running `Metabolomics2Network`. Per
# default, it is `chebi class mapping`.
# 4. `distance`. Distance in the ontology between the original metabolite and the
# metabolite to which it was mapped.
# 5. `chebi`. ChEBI id of the metabolite to which the original metabolite was
# mapped.


#' @name .mapMetToGSMN_metabolomics2network
#'
#' @aliases mapMetToGSMN_metabolomics2network
#'
#' @title Map the experimental nodes to GSMN nodes using Metabolomics2Network Python tool.
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
.mapMetToGSMN_metabolomics2network <- function(inputData, resFile = "Res_Met2Net_MappedMet.txt") {

    ## checks
    if(inputData$met2NetDir==""){print("met2NetDir is empty");return()}
    if(inputData$metF==""){print("metF is empty");return()}
    if(inputData$resPath==""){print("resPath is empty");return()}
    if(resFile==""){print("resFile is empty");return()}
    if(inputData$configF==""){print("configF is empty");return()}

    pathToMappings <- paste0(inputData$resPath, "GSMNMappings/")
    dir.create(pathToMappings)

    # generate command line to execute metabolites2Network
    com <-
      paste0(
        "python3 ",
        inputData$met2NetDir,
        "metabolomics2network.py",  # Python package file
        " tsv ",                    # file_type
        inputData$idenMetF,         # metabolomics_path
        " ",
        inputData$metF,            # network_metabolites_path
        " ",
        pathToMappings,          # output_path
        resFile,
        " ",
        inputData$configF,          # conf_file_path
        " 1,2"                        # mapping_types, 1: exact multimapping, 2: chebi class mapping
        )

    # execute code
    system(com)

    # process the mappings to clean the results
    processMappings(identMetF = inputData$idenMetF, pathToMappings, resFile)

    return()

}



#' @name .mapMetToGSMN_inchikey
#'
#' @title Map features to metabolites if the InChIKey is the same.
#'
#' @description
#' Function to map the identified experimental nodes to the
#' corresponding GSMN nodes, by InChIKey comparison.
#'
#' @param inputData
#' `list`, list returned by the `loadInputData` function
#'
#' @return
#'
#' @author Sarah Scharfenberg
#'
.mapMetToGSMN_inchikey <- function(inputData){

  # get metadata from QFeatures object
  metadata <- rowData(inputData$peakList[[1]])

  ## metadata$database_identifier
  ## annotierte chebi ID

  if(sum(!is.na(metadata$inchi))==0){ print("No inchis annotated."); return()}


}
