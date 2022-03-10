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
#'  @param assay_name
#'  `Character`, define which assay needs to be extracted from QFeature input
#'  e.g. "pos".
#'
#' @param ...
#'  Insert here parameter from `structural` function from `MetNet`package.
#'  (https://github.com/MetClassNet/MetNet)
#'  `transformation` is a  `data.frame`, containing the columns `"group"`,
#'  and `"mass"` that will be used for detection of transformation of
#'  (functional) groups
#'  Parameter `ppm`is `numeric`, mass accuracy of m/z features in
#'  parts per million (ppm)
#'
#' @import
#' `MetNet`
#' `QFeatures`
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
qfeat_structural <- function(x, assay_name = "features", ...) {

  feat_int <- as.data.frame(assay(x[[assay_name]]))
  feat_names <- as.data.frame(rowData(x[[assay_name]]))

  feat_names <- feat_names[,c("metabolite_identification", "mz", "rtime")]
  colnames(feat_names) <- c("manualAnnotation", "mz", "RT")

  feat <- merge(feat_int, feat_names, by = "row.names", sort = FALSE, all=TRUE)
  row.names(feat) <- feat$"Row.names"
  feat <- subset(feat, select = - c(Row.names))

  MetNet::structural(x = feat, ...)

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
#'  @param assay_name
#'  `Character`, define which assay needs to be extracted from QFeature input
#'  e.g. "pos".
#'
#'  @param ...
#'  Insert here parameter from `statistical` function from `MetNet`package.
#'  (https://github.com/MetClassNet/MetNet)
#'  `model` is a `character` vector containing the methods that will be used
#'  (`"lasso"`, `"randomForest"`, `"clr"`, `"aracne"`, `"pearson"`,
#'  `"pearson_partial"`, `"pearson_semipartial"`, `"spearman"`,
#'  `"spearman_partial"`, `"spearman_semipartial"`, `"bayes"`)
#'  `data.frame`, containing the columns `"group"`.
#'  `p`is `logical`, by default is set to FALSE. `
#'
#' @import
#' `MetNet`
#' `QFeatures`
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
qfeat_statistical <- function(x, assay_name = "features",
                              na.omit = FALSE, ...) {

   feat_int <- as.matrix(assay(x[[assay_name]]))

  if (na.omit == TRUE) {
    feat_int <- feat_int %>% na.omit() %>% as.matrix()
  }


  MetNet::statistical(feat_int, ...)

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
#' @param assay_name
#'  `Character`, define which assay needs to be extracted from QFeature input
#'  e.g. "pos".
#'
#' @param ...
#'  Insert here parameter from `homol.search` function from `nontarget`package.
#'  (https://cran.r-project.org/web/packages/nontarget/index.html)
#'  Define parameters needed at `homol.search` as `elements`, `use_C`,
#'  `minmx`, `maxmz`, `minrt`, `maxrt`, `ppm`, `mztol`, `rttol`, `minlenght`,
#'  `mzfilter`, `spar`, `R2`.
#'  For further information see Help page of `?homol.search`.
#'
#' @import
#' `MetNet`
#' `QFeatures`
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
qfeat_homol <- function(x, assay_name = "features", plot = FALSE, ...) {

  feat_int <- as.data.frame(assay(x[[assay_name]]))
  feat_names <- as.data.frame(rowData(x[[assay_name]]))

  ## Extract feature definitions
  featid <- feat_names$database_identifier
  rt <- feat_names$rtime
  mz <- feat_names$mz

  data(isotopes)

  ## Use peaklist
  peaklist <- data.frame(mass=mz,
                         intensity=feat_int[,1],
                         rt=rt)

  homol <- nontarget::homol.search(peaklist,
                                   isotopes,
                                   ...)

  #(4.2) Plot results
  if(plot) {

    nontarget::plothomol(homol,
                         xlim=FALSE,
                         ylim=FALSE,
                         plotlegend=TRUE)

  }

  homol

}
