#' @example 
#' library(QFeatures)
#' library(tidyverse)
#' data_df <- read_tsv("/Users/Liesa4/Library/Mobile Documents/com~apple~CloudDocs/Promotion/R/Projects/MTBLS291/m_mtbls291_POS_mass_spectrometry_v2_maf.tsv")
#' data_qF <- readQFeatures(data_df, ecol = 22:ncol(data_df), name = "pos")
#' transformations <- read.csv('/Users/Liesa4/Library/Mobile Documents/com~apple~CloudDocs/Promotion/R/Projects/MTBLS291/Lipid_table.csv') 
#' %>% select(group, formula, mass) %>% as.data.frame()
#'
#' @param x `QFeature`
#'
qfeat_structural <- function(x, transformation, ...) {
  feat_int <- assay(x[["pos"]]) %>% as.data.frame()
  feat_names <- rowData(x[["pos"]]) %>% as.data.frame() %>% 
    dplyr::rename(mz = mass_to_charge) %>% 
    dplyr::rename(manualAnnotation = `metabolite_identification` ) %>%
    dplyr::rename(RT = retention_time) %>%
    dplyr::select(manualAnnotation, mz, RT)
  
  feat <- merge(feat_int, feat_names, by = "row.names", sort = F, all=T ) %>% 
    select(- `Row.names`)
  
  mass_diff <- MetNet::structural(x = feat, transformation = transformations, ppm = 10)
  
  return(mass_diff)
}

#'
#'
#' @param x `QFeature`
#'
qfeat_statistical <- function(x, model, p = F, ...) {
  feat_int <- assay(x[["pos"]]) %>% as.matrix() 
  
  if (p == T){
    if (model == "spearman" || model == "pearson" ) {
      corr <- MetNet::statistical(feat_int, model = model, p = TRUE)
      
    }
    else if (model != "spearman" || model != "pearson"){
      stop("'model' not implemented in statistical(p = TRUE)")
      
    }
  }
  else if (p == FALSE){
    corr <- MetNet::statistical(feat_int, model = model, p = F)
  }
  
  
  return(corr)
}

#'
#'
#' @param x `QFeature`
#'
qfeat_homol <- function(x, ...) {

}
