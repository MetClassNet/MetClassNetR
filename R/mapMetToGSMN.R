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
.mapMetToGSMN_inchikey <- function(inputData, gsmnMetaFile=NULL){

  ## get inchikeys from QFeatures object
  metadata <- rowData(inputData$peakList[[1]])

  if(sum(!is.na(metadata$inchi))==0){
    print(paste0("No InChIKeys annotated to the experiment. Please add InChIKeys to your metadata or",
          " try annotate_chemical_properties() before calling again."))
    return()
  }else{
    data_IK <- cbind(metadata$id,metadata$inchi)
  }

  ## get inchis from gsmn
  if(is.null(gsmnMetaFile)){

    gsmn_data <- read.table(inputData$metF, sep="\t", header=TRUE, stringsAsFactors = FALSE,quote="", comment.char = "")
    if(any(grepl("inchi",names(gsmn_data), ignore.case = TRUE))){
      print(paste0("No InChIKeys annotated to the GSMN. Please add InChIKeys to your metadata or",
                   "provide the information in a separate file. To retrieve chemical properties based on ChEBI Ids",
                   " use getChemicalPropertiesFromChebi() to create such a file."))
      return()
    }else{
      ikCol <- which(grepl("inchi",names(gsmn_data), ignore.case = TRUE))
      gsmn_IK <- cbind(gsmn_data$ID,gsmn_data[,ikCol])
    }
  }else{
    gsmn_data <- read.table(gsmnMetaFile, sep="\t", header=TRUE, stringsAsFactors = FALSE,quote="", comment.char = "")
    gsmn_IK <- cbind(gsmn_data$ID,gsmn_data[,ikCol])

  }


  ## compare per inchikey


  ## write results in mapping file
 # <- c("metabolite.name","mapped.on.id","mapping.types","distance","chebi","originalChebi"
#                   Cluster_0497	M_eicostet	chebi class mapping	-2	['82835']	CHEBI:132539



}


#' @name getChemicalPropertiesFromChebi
#'
#' @title Function to retrieve chemical properties such as InChIKey, Smiles
#' from a GSMN file containing ChEBI Ids.
#'
#' @export
getChemicalPropertiesFromChebi <- function(fileIn, fileOut=NULL){

  ## get inchikeys from igraph object of the GSMN
  gsmn_data <- read.table(source, sep='\t', header=TRUE, stringsAsFactors = FALSE,quote="", comment.char = "")
  n <- length(unique(gsmn_data$Chebi))

  gsmn_metadata <- sapply(unique(gsmn_data$Chebi), function(chebi_id){
    prop <- webchem::chebi_comp_entity(chebi_id)[[1]]$properties
    return(c(chebi_id,prop$inchi,prop$inchikey,prop$smiles))
    })
  gsmn_metadata <- t(gsmn_metadata)
  colnames(gsmn_metadata) <- c("chebiid","inchi","inchikey","smiles")

  res <- cbind(gsmn_data,gsmn_metadata[gsmn_data$Chebi,2:4])

  if(!is.null(fileOut)){
     write.table(res,
                file=fileOut,
                row.names = FALSE,
                col.names = TRUE,
                sep='\t')
  }
  return(res)
}

#' @name annotateChemicalPropertiesToExperiment
#'
#' @title Function to annotate InputData with chemical properties.
#'
#' @description This function retrieves chemical properties such as SMILES, InChI and InChIKey
#' from ChEBI using the package webchem. Notice that there are ChEBI entries which have SMILES but
#' no InChI/InChIKey annotated because of the possibility to represent variable ends with an *
#' in SMILES but not in InChI/InChIKey.
#'
#' @param InputData QFeatures experiment
#'
#' @return QFeatures experiment similar to input and with chemical properties added to peakList RowData
#'
#' @author Sarah Scharfenberg
#'
#' @export

annotateChemicalPropertiesToExperiment <- function(inputData){

  metadata <- rowData(inputData$peakList[[1]])

  if(sum(!is.na(metadata$database_identifier))==0){
    print("No database_identifier found to retrieve InChIKey from.")
    return()
  }else{
    print("Trying to use database_identifier to retrieve InChIKeys from ChEBI.")

   any_database_id <- metadata$database_identifier[metadata$database_identifier!=""][1]
   if(!grepl("CHEBI:", any_database_id)){
     print("database_identifyer does not contain ChEBI IDs")
     return()
   }else{
    print("Found ChEBI IDs in database_identifier.")

    all_entities <- webchem::chebi_comp_entity(metadata$database_identifier[metadata$database_identifier!=""])

    all_inchikey <- sapply(metadata$database_identifier, function(id){ return(all_entities[[id]]$properties$inchikey)})
    ## if there is no chebi id, the call will return NA, if there is a chebi id but no inchi annotated the call will return null
    ## both are reported as empty strings in the result
    all_inchikey <- lapply(all_inchikey,function(inchikey){if(is.null(inchikey)){return("")}else{return(inchikey)}})
    all_inchikey <- lapply(all_inchikey,function(inchikey){if(is.na(inchikey)){return("")}else{return(inchikey)}})
    all_inchikey <- unlist(all_inchikey)

    all_inchi <- sapply(metadata$database_identifier, function(id){ return(all_entities[[id]]$properties$inchi)})
    all_inchi <- lapply(all_inchi,function(inchi){if(is.null(inchi)){return("")}else{return(inchi)}})
    all_inchi <- lapply(all_inchi,function(inchi){if(is.na(inchi)){return("")}else{return(inchi)}})
    all_inchi <- unlist(all_inchi)

    all_smiles <- sapply(metadata$database_identifier, function(id){ return(all_entities[[id]]$properties$smiles)})
    all_smiles <- lapply(all_smiles,function(smiles){if(is.null(smiles)){return("")}else{return(smiles)}})
    all_smiles <- lapply(all_smiles,function(smiles){if(is.na(smiles)){return("")}else{return(smiles)}})
    all_smiles <- unlist(all_smiles)

    dF <- DataFrame(inchi = all_inchi,
                    inchikey = all_inchikey,
                    smiles = all_smiles)

    rowData(inputData$peakList) <- List(features = dF)
       }
  }

    return(inputData)
}
