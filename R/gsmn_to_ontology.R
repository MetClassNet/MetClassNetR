#'@title Ontology classification for a whole network.
#'@description This function will annotate all compounds of the given GSMN
#'with a class of the chosen ontology.
#'
#'@usage gsmn_to_ontology(file=NULL, ontology = c("chemont","chebi"),
#'output=NULL, delay=0)
#'
#'@param file SBTab file of a genome scale metabolic network (GSMN) holding a unique
#'identifier (column "ID") and an inchikey (column "Identifiers.inchikey")
#'@param ontology character that defines which ontology is used for the mapping.
#'Can be one of c("chemont", "chebi")
#'@param output optional, a filename for a two column output table
#'@param delay numeric time delay bewtween two requests in classyfireR. See details.
#'
#'@details
#'For a given SMILES or InChIKey the package classyfireR will report the classes
#'a compound is classified into in ChemOnt. A lot of compounds are already classified.
#'Still there remain some unclassified compounds that appear in this script as NULL
#'or empty in the result list. These queries can be filtered and manually requested
#'for a new classification via the Webservice available at http://classyfire.wishartlab.com/.
#'
#' When requesting a high number of queries within classyfireR the following error
#' message appears: "Error in classyfireR::get_classification(ik) : Request rate limit
#' exceeded!" That refers to a limit of requests per time. To address that set
#' delay=10 meaning that the script will wait 10 seconds after each request.
#' The result will take much longer but will finish in the end.
#'
#'@note
#'
#'@return A csv file giving the most detailed classification per compound.
#'
#'@author Sarah Scharfenberg
#'
#'@examples
#'file <- system.file("extdata/Example_Coverage/iHY3410_Compound-SBtab.tsv", package = "MetClassNetR")
#'gsmn_to_ontology(file,"chemont", delay=10)
#'
#'@export
#'
gsmn_to_ontology <- function(file=NULL,
                             ontology = "chemont",
                             output = NULL,
                             delay=0){

  # Parameter checks
  if(is.null(file)){
    print("No file provided.")
    return(NULL)
  }
  if(!file.exists(file)){
    print(paste0("File ",file, " not found."))
    return(NULL)
  }
  if(!is.element(ontology,c("chemont","chebi"))){
    print("Please choose one of chebi or chemont as ontology.")
    print("To include other ontologies feel free to generate a proper pull request.")
    return(NULL)
  }
  ##


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

  network_tibble <- dplyr::as_tibble(network[c("ID","Identifiers.inchikey")])

  if(ontology == "chemont"){
    network_tibble <- network_tibble %>%
      get_classification_chemont(delay=delay)
  }
  if(ontology == "chebi"){
    network_tibble <- network_tibble %>%
      get_classification_chebi()
  }

  if(!is.null(output)){
    network_tibble %>%
      select("ID","class") %>%
      readr::write_csv(file=output)
  }

  return(network_tibble)
}
