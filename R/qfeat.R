#' @name qfeat_structural
#' 
#' @aliases qfeat_structural
#' 
#' @title Create m/z difference adjacency matrix from QFeatures input
#' 
#' @description 
#' The function `qfeat_structural` uses the input from `QFeatures` and creates an 
#' adjacency matrix based on m/z (molecular weight) difference from the `MetNet` package. 
#' 
#' @param x 
#' `QFeatures` file
#' 
#' @param transformation  
#' `data.frame`, containing the columns `"group"`,
#' and `"mass"` that will be used for detection of transformation of
#' (functional) groups (https://github.com/MetClassNet/MetNet)
#' 
#' @param ppm
#' `numeric`, mass accuracy of m/z features in parts per million (ppm)
#' @import 
#' `MetNet` 
#' `QFeatures`
#' `dplyr`
#' 
#' @details 
#' `qfeat_structural` extracts required information from a `QFeatures` input 
#'  and builds data.frames containing intensity informations, m/z values, 
#'  RT values, and annotations of all samples. 
#'  Then the function `structural` from the `MetNet` package is applied
#'  to generate a m/z difference adjacency matrix. 
#'  
#'  @return 
#'  `list` containing three matrices. The first entry stores the `numeric`
#' `matrix` with edges inferred from mass differences. The second entry
#' stores the `character` `matrix` with the type (corresponding to the
#' `"group"` column in `transformation`) is stored. The third matrix stores
#' another `character` `matrix` with the corresponding mass difference (`"mass"`
#' column in `transformation`).
#' 
#' @author Liesa Salzer,  \email{liesa.salzer@@helmholtz-muenchen.de}
#' 
#' @examples 
#' ####### To be added
#' 
#' @export
qfeat_structural <- function(x, transformation, ppm = 10, ...) {
  feat_int <- assay(x[["pos"]]) %>% as.data.frame()
  feat_names <- rowData(x[["pos"]]) %>% as.data.frame() %>% 
    dplyr::rename(mz = mass_to_charge) %>% 
    dplyr::rename(manualAnnotation = `metabolite_identification` ) %>%
    dplyr::rename(RT = retention_time) %>%
    dplyr::select(manualAnnotation, mz, RT)
  
  feat <- merge(feat_int, feat_names, by = "row.names", sort = F, all=T ) %>% 
    select(- `Row.names`)
  
  mass_diff <- MetNet::structural(x = feat, transformation = transformations, ppm = ppm)
  
  return(mass_diff)
}




#' @name qfeat_statistical
#' 
#' @aliases qfeat_statistical
#' 
#' @title Create a list of correlation matrices from QFeatures input
#' 
#' @description 
#' The function `qfeat_statistical` uses the input from `QFeatures` and creates an 
#' adjacency matrix based on statistical methods using the `MetNet` package. 
#' The function includes functionality to calculate adjacency matrices based on
#' LASSO (L1 norm)-regression, random forests, context likelihood of
#' relatedness (CLR), the algorithm for the reconstruction of accurate
#' cellular networks (ARACNE), Pearson correlation (also partial and
#' semipartial), Spearman correlation (also partial and semipartial)
#' and score-based structure learning (Bayes). The function returns a
#' list of adjacency matrices that are defined by `model`.
#' Additionally, for pearson and/or spearman correlation also the negative
#' correlation values and the corresponding p-Value is calculated and listed,
#' when `p` is set to TRUE.
#' 
#' @param x 
#' `QFeatures` file
#' 
#' @param model  
#' `character` vector containing the methods that will be used
#' (`"lasso"`, `"randomForest"`, `"clr"`, `"aracne"`, `"pearson"`,
#' `"pearson_partial"`, `"pearson_semipartial"`, `"spearman"`,
#' `"spearman_partial"`, `"spearman_semipartial"`, `"bayes"`)
#' `data.frame`, containing the columns `"group"`
#' (https://github.com/MetClassNet/MetNet).
#' 
#' @param p
#' `logical`, by default is set to FALSE. 
#' 
#' @import 
#' `MetNet` 
#' `QFeatures`
#' `dplyr`
#' 
#' @details 
#' `qfeat_statistical` extracts required information from a `QFeatures` input 
#'  and builds data.frames containing intensity informations of all samples. 
#'  Then the function `statistical` from the `MetNet` package is applied
#'  to calculate adjacency matrices based on
#'  LASSO (L1 norm)-regression, random forests, context likelihood of
#'  relatedness (CLR), the algorithm for the reconstruction of accurate
#'  cellular networks (ARACNE), Pearson correlation (also partial and
#'  semipartial), Spearman correlation (also partial and semipartial)
#'  and Constraint-based structure learning (Bayes).
#'  The default of `p` is FALSE. Then all types can be selected in `model`. 
#'  The positive correlation value will be displayed in the correlation matrix. 
#'  If `p` is set to TRUE, only "pearson" and/or "spearman" correlation may be 
#'  selected. As output positive and negative correlation values will be displayed. 
#'  Morover their corresponding p-values will be added to the list. 
#'  
#'  @return 
#'  `list` containing the respective adjacency matrices specified by
#' `model`.  It p` is TRUE, also the corresponding p-values of Spearman and/or
#'  Pearson Correlation will be added to the `list`.
#' 
#' @author Liesa Salzer,  \email{liesa.salzer@@helmholtz-muenchen.de}
#' 
#' @examples 
#' ####### To be added
#' 
#' @export
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

#' @name qfeat_homol
#' 
#' @aliases qfeat_homol
#' 
#' @title Homologue series extraction from QFeatures input
#' 
#' @description 
#' The function `qfeat_homol` uses the input from `QFeatures` and uses dynamic 
#' programming algorithms from `homol.search` from the `nontarget` package 
#' (https://cran.r-project.org/web/packages/nontarget/index.html) for
#' unsupervides detection of homologue series in LC-(HR)MS data.
#' 
#' @param x 
#' `QFeatures` file
#' 
#' @param elements 
#' FALSE or chemical elements in the changing units of the homologue 
#' series, e.g. c("C","H") for alkane chains. Used to restrict search.
#' By default set to c("C","H","O").
#' 
#' @param use_C
#' For elements: take element ratio to C-atoms into account? Used to restrict search.
#' By default set to TRUE.
#' 
#'  @param minmz 
#'  Defines the lower limit of the m/z window to search homologue series peaks, 
#'  relative to the m/z of the one peak to search from. Absolute m/z value [u].
#'  
#'  @param maxmz
#'  Defines the upper limit of the m/z window to search homologue series peaks, 
#'  relative to the m/z of the one peak to search from. Absolute m/z value [u].
#'  
#'  @param minrt
#'  Defines the lower limit of the retention time (RT) window to look for 
#'  other ho- mologue peaks, relative to the RT of the one peak to search from, 
#'  i.e., RT+minrt. For decreasing RT with increasing HS mass, use negative 
#'  values of minrt.
#'  
#'  @param maxrt
#'  Defines the upper limit of the RT window to look for other homologue peaks, 
#'  relative to the RT of the one peak to search from, i.e., RT+maxrt. See `minrt`.
#'  
#'  @param ppm
#'  Should mztol be set in ppm (TRUE, default) or in absolute m/z [u] (FALSE)?
#'  
#'  @param mztol 
#'  m/z tolerance setting: +/- value by which the m/z of a peak may vary from 
#'  its expected value. If parameter ppm=TRUE (see below) given in ppm, 
#'  otherwise, if ppm=FALSE, in absolute m/z [u]. Default is 5.
#'  
#'  @param rttol
#'  Retention time (RT) tolerance by which the RT between two adjacent pairs of 
#'  a homologue series is allowed to differ. Units as given in column 3 of 
#'  peaklist argument, e.g. [min]. Default is 5. 
#'  
#'  @param minlength
#'  Minimum number of peaks in a homologue series. Default is 4.
#'  
#'  @param mzfilter
#'  Vector of numerics to filter for homologue series with specific m/z 
#'  differences of their repeating units, given the tolerances in mztol. 
#'  Mind charge z! Default is FALSE. 
#'  
#'  @param spar
#'  Smoothing parameter, typically (but not necessarily) in (0,1).
#'  Default 0.45. 
#'  
#'  @param R2
#'  FALSE or 0<numeric<=1. Coefficient of determination for cubic smoothing 
#'  spline fits of m/z versus retention time; homologue series with lower R2 are 
#'  rejected. See smooth.spline. Default 0.98. 

#'  
#' @import 
#' `MetNet` 
#' `QFeatures`
#' `dplyr`
#' 
#' @details 
#' `qfeat_structural` extracts required information from a `QFeatures` input 
#'  and builds `list` containing results from `homol.search` function of 
#'  `nontarget` package. 
#'  
#'  @seealso `nontarget`, `homol.search`
#'  https://cran.r-project.org/web/packages/nontarget/index.html
#'  
#'  @return 
#'  `list` of type homol with 6 entries. First entry contains dataframe 
#'  with peaks (mass,intensity,rt,peak ID) and their homologue 
#'  series relations (to ID,m/z increment,RT increment) within different 
#'  homologue series (HS IDs,series level). Last column HS cluster states 
#'  HS clusters into which a peak was assigned via its HS. Second entry 
#'  contains the parameters used.
#'  The third list entry is a Dataframe listing all peaks (peak IDs) per 
#'  homologue series (HS IDs), the underlying mean m/z & RT increments 
#'  (m/z increments, RT increments) and the minimum and maximum RT changes 
#'  between individual peaks of the series. The 4th entry stores used 
#'  m/z restrictions. List of peak IDs per level in the individual series
#'  is stores in the 5th entry. The 6th entry contains List with 
#'  superjacent HS IDs per group - for setdeb=c(3,...).
#'  
#'  
#' 
#' @author Liesa Salzer,  \email{liesa.salzer@@helmholtz-muenchen.de}
#' 
#' @examples 
#' ####### To be added
#' 
#' @export
qfeat_homol <- function(x, 
                        elements=c("C","H","O"), use_C=TRUE,
                        minmz, 	maxmz,
                        minrt,  maxrt,
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


