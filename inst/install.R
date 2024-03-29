install.packages("devtools")
library("devtools")

##
## CRAN
##
install.packages(c("nontarget", "Hmisc",
	"dplyr", "igraph", "networkD3", "tibble", "tidyverse", "visNetwork",
	"rmarkdown", "knitr", "ggraph", "rcdk", "VennDiagram", "ontologyIndex"))

##
## BioC
##
BiocManager::install(c("HDF5Array", "BiocStyle",
	"MsCoreUtils", "Spectra", "QFeatures", #"MsBackendMgf",
	"MSnbase"))

##
## Special stuff
##
install_github("aberHRML/classyfireR")
install_github("tnaake/MetNet")
install_github("rformassspectrometry/MsBackendMgf")
install_github("sneumann/MsBackendMsp")
install_github("CDK-R/rinchi")


## Should be in BioC ?!
#install_github("rformassspectrometry/QFeatures")
#install_github("rformassspectrometry/MsBackendMgf")
