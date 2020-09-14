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
qfeat_homol <- function(x, 
                        elements=c("C","H","O"), use_C=TRUE,
                        minmz=5, 	maxmz=50,
                        minrt=5,  maxrt=60,
                        ppm=TRUE,
                        mztol=5,  rttol=5,
                        minlength=4,
                        mzfilter=FALSE,
                        spar=.45, 	R2=.98, ...) {
  
  feat_int <- assay(x[["pos"]]) %>% as.data.frame()
  feat_names <- rowData(x[["pos"]]) %>% as.data.frame()
  
  ## Extract feature definitions
  featid <- feat_names$database_identifier
  rt <- feat_names$retention_time
  mz <- feat_names$mass_to_charge
  
  # (0.2) list of isotopes - package enviPat ############
  data(isotopes)
  
  ## Use peaklist
  peaklist <- data.frame(mass=mz, intensity=feat_int[,1], rt=rt)
  
  
  homol <- nontarget::homol.search(peaklist,
                                   isotopes,		elements=elements, use_C=use_C,
                                   minmz=minmz, 	maxmz=maxmz,
                                   minrt=minrt,  maxrt=maxrt,
                                   ppm=ppm,
                                   mztol=mztol,  rttol=rttol,
                                   minlength=minlength,
                                   mzfilter=mzfilter,
                                   spar=spar, 	R2=R2,
                                   plotit=FALSE)
  
  #(4.2) Plot results 
  nontarget::plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
  
  return(homol)
}


