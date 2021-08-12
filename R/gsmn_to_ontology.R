#' @name gsmn_to_ontology
#'
#' @title Ontology classification for a whole network
#'
#' @description This function will annotate all compounds of the given GSMN
#' with a class of the chosen ontology.
#'
#' @usage gsmn_to_ontology(file=NULL, inchikeys=NULL, smiles=NULL, inchis=NULL,
#' ontology = c("chemont","chebi"), output=NULL, delay=0, conn=NULL)
#'
#' @param file
#'  SBTab file of a genome scale metabolic network (GSMN) holding a unique
#'  identifier (column "ID") and an inchikey (column "Identifiers.inchikey"). Only one
#'  of the parameters file, inchikeys, smiles or inchis should be provided.
#'
#' @param inchikeys
#'  a vector of InChIKeys that will be be classified. Only one of the
#'  parameters file, inchikeys, smiles or inchis should be provided.
#'
#' @param smiles
#'  a vector of SMILES that will be classified. Only one of the
#'  parameters file, inchikeys, smiles or inchis should be provided.
#'
#' @param inchis
#'  a vector of InChIs that will be converted. Only one of the
#'  parameters file, inchikeys, smiles or inchis should be provided.
#'
#' @param ontology
#'  character that defines which ontology is used for the mapping.
#'  Can be one of c("chemont", "chebi")
#'
#' @param output
#'  optional, a filename for a two column output table
#'
#' @param delay
#'  numeric time delay between two requests in classyfireR. See details.
#'
#' @details
#' For a given SMILES or InChIKey the package classyfireR will report the classes
#' a compound is classified into in ChemOnt. A lot of compounds are already classified.
#' Still there remain some unclassified compounds that appear in this script as NULL
#' or empty in the result list. These queries can be filtered and manually requested
#' for a new classification via the Webservice available at http://classyfire.wishartlab.com/.
#'
#' When requesting a high number of queries within classyfireR the following error
#' message appears: "Error in classyfireR::get_classification(ik) : Request rate limit
#' exceeded!" That refers to a limit of requests per time. To address that set
#' delay=10 meaning that the script will wait 10 seconds after each request.
#' The result will take much longer but will finish in the end.
#'
#' @note
#'
#' @return A csv file giving the most detailed classification per compound.
#'
#' @author Sarah Scharfenberg
#' @export
#'
#' @examples
#' \dontrun{
#' file <- system.file("extdata/Example_Coverage/iHY3410_Compound-SBtab.tsv", package = "MetClassNetR")
#' gsmn_to_ontology(file,"chemont", delay=10)
#' }
#'
gsmn_to_ontology <- function(file=NULL,
                             inchikeys=NULL,
                             smiles = NULL,
                             inchis = NULL,
                             ontology = "chemont",
                             output = NULL,
                             delay=0,
                             conn=NULL){

  ## Check Input Parameter
  if(sum(c(!is.null(file), !is.null(inchikeys), !is.null(smiles), !is.null(inchis)))!=1){
    print("Wrong input provided. Please provide only one of file, inchikey, smiles or inchis.")
    return(NULL)
  }
  if(!is.null(file)){
    if(!file.exists(file)){
    print(paste0("File ",file, " not found."))
    return(NULL)
  }}
  if(!is.element(ontology,c("chemont","chebi"))){
    print("Please choose one of chebi or chemont as ontology.")
    print("To include other ontologies feel free to generate a proper pull request.")
    return(NULL)
  }


  ## Prepare Input from file
  if(!is.null(file)){
    ## Read SBTab to tibble
    network <- read.csv(file,
                        header = TRUE,
                        skip=1,
                        stringsAsFactors = FALSE,
                        sep="\t",
                        na.strings = "")

    colnames(network) <- gsub(pattern = "X.",
                              replacement = "",
                              x = colnames(network))

    network_tibble <- na.omit(dplyr::as_tibble(network[c("ID","Identifiers.inchikey")]))
  }

  ## Prepare Input from smiles vector
  if(!is.null(smiles)){

    network_tibble <- new_tibble(data.frame(ID=smiles)) %>% na.omit()

    netinchiks <- data.frame(sapply(network_tibble$ID,function(smi){
      if(!is.na(smi)){
        smim <- rcdk::parse.smiles(smi)
        if(!is.null(unlist(smim))){
          res <- rinchi::get.inchi.key(smi)
          return(res)
        }else{return(NA)}
      }else{return(NA)}
    }))

    print("The following smiles could not be converted:")
    print(paste(network_tibble$ID[is.na(netinchiks)], sep="\n"))

    network_tibble <- network_tibble %>%
      add_column(Identifiers.inchikey=as.character(unlist(netinchiks))) %>%
      na.omit()
  }

  ## Prepare Input from inchi vector
  if(!is.null(inchis)){

    network_tibble <- new_tibble(data.frame(ID=inchis)) %>% na.omit()

    netinchiks <- data.frame(sapply(network_tibble$ID,function(inchi){
      if(!is.na(inchi)){
        inchim <- rcdk::parse.inchi(inchi)
        if(!is.null(unlist(inchim))){
          res <- rinchi::get.inchi.key(inchim)
          return(res)
        }else{return(NA)}
      }else{return(NA)}
    }))

    network_tibble <- network_tibble %>%
      add_column(Identifiers.inchikey=as.character(unlist(netinchiks))) %>%
      na.omit()

  }


  ## Prepare Input from inchikey vector
  if(!is.null(inchikeys)){
    network_tibble <- new_tibble(data.frame(ID=inchikeys, Identifiers.inchikey=inchikeys)) %>%
                        na.omit()
  }

  ## Classification
  if(ontology == "chemont"){
    network_tibble <- network_tibble %>%
      get_classification_chemont(delay=delay, conn=conn)
  }
  if(ontology == "chebi"){
    network_tibble <- network_tibble %>%
      get_classification_chebi()
  }


  ## Write and Report results
  if(!is.null(output)){
    network_tibble %>%
      select("ID","class") %>%
      readr::write_csv(file=output)
  }

  return(network_tibble)
}
