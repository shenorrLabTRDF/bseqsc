# Functions to perform the deconvolution analysis
# 
# Author: Renaud Gaujoux
###############################################################################

#' Scale Single Cell Data According to Average Cell Type Profiles
#' 
#' @param x matrix of raw counts
#' @param clusters vector giving the cell type of each cell (i.e. column) in the data `x`.
#' If `x` is an `ExpressionSet` object, it can also be the name (as a string) of the 
#' phenotypic variable to use.
#' @param samples vector giving the sample of origin of each cell (i.e. column) in the data `x`.
#' If `x` is an `ExpressionSet` object, it can also be the name (as a string) of the 
#' phenotypic variable to use.
#' 
cpm_cell_type <- function(x, clusters, samples){
	
  # extract variables
  clusters <- pVar(x, clusters)
  samples <- pVar(x, samples)
  
	# compute total count
	tot <- colSums(exprs(x))
	# compute CPM
  cpm <- sweep(exprs(x), 2L, tot, '/')
	# re-scale with within cell type average
	sc <- as.character(paste(clusters, samples))
	tot_map <- sapply(split(tot, sc), mean)
	cpm <- sweep(cpm, 2L, tot_map[sc], '*')
	
  # result
  res <- x
  exprs(res) <- cpm
  res
}

#' Building a Basis Matrix of Reference Gene Expression Profiles
#' 
#' 
#' @inheritParams cpm_cell_type
#' @param ct.scale logical that indicates if the single cell expresson profiles
#' should be rescaled according to their cell type average total count.
#' @param markers list of cell type marker genes.
#' The type of gene identifiers must be the same as the ones used as feature/row names in `x`.
#' @param limit logical that indicates if the returned basis matrix should only contain
#' cell types listed in `markers`.
#' 
#' @return a `matrix` object.
#' 
#' @export
bseqsc_basis <- function(x, markers, clusters, samples, ct.scale = TRUE, limit = TRUE){
  
  # scale data before subsetting to markers
  if( ct.scale ){
    x <- cpm_cell_type(x, clusters = clusters, samples = samples)
  }
  
  # compute averages on markers
  ids <- intersect(unlist(markers), rownames(x))
  x <- x[ids, , drop = FALSE]
  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))
  # within each cell type
  res <- sapply(unique(clusters), function(ct){
      # mean of sample means
      rowMeans(sapply(unique(samples), function(sid){
        y <- exprs(x)[, clusters %in% ct & samples %in% sid, drop = FALSE]
        rowMeans(y)
      }), na.rm = TRUE)
    })

  if( limit ) res <- res[, names(markers)]
  res
}

#' Estimating Cell Type Proportions From Single Cell Data
#' 
#' @param x Expression data from heterogenous samples
#' @param reference Bisis matrix of reference cell type expression profiles.
#' If `NULL` then a built-in reference for pancreatic islet cell sub-population
#' is used (data [PancreasIslet]).
#' @param log logical that indicates if the data `x` is in in log-scale.
#' If `NULL`, then log scale is inferred by [xbioc::is_logscale].
#' @param ... other arguments passed to `CIBERSORT`.
#' @param verbose logical that toggles log messages.
#' 
#' @return a list with elements:
#'   * coef: matrix of estimated proportions (cell type x samples)
#'   * stats: statistics computed by `CIBERSORT`
#' 
#' @export
#' @importFrom e1071 svm
#' @importFrom parallel mclapply
#' @importFrom preprocessCore normalize.quantiles
bseqsc_proportions <- function(x, reference = NULL, log = NULL, ..., verbose = TRUE){
	
  message <- if( verbose ) message else function(...) invisible(NULL)
  
  # load pancreas basis matrix if necessary
  if( is.null(reference) ){
    message("* Using pancreatic islet reference basis matrix: ", appendLF = FALSE)
    reference <- ldata('PancreasIslet')
  }
  
	# load CIBERSORT
  CIBERSORT <- bseqsc_config(error = TRUE)
	y <- x
  if( is(y, 'ExpressionSet') ) y <- exprs(y)
  
  # use common features only
  ids <- intersect(rownames(y), rownames(reference))
  message(sprintf("* Data features: %s", str_out(rownames(y), total = TRUE)))
  message(sprintf("* Basis features: %s", str_out(rownames(reference), total = TRUE)))
  message(sprintf("* Common features: %s", str_out(ids, total = TRUE)))
  islog <- log %||% is_logscale(y)  
  if( islog ){
    message("* Converting to linear scale")
    y <- expb(y, 2)
  }
  
  # setup run directory and files
  tdir <- tempfile('CIBERSORT')
  dir.create(tdir)
  owd <- setwd( tdir )
  on.exit({
        unlink(tdir, recursive = TRUE)
        setwd(owd)
      }) 
  message("* Writing input files ... ", appendLF = FALSE)
  write.table(reference, file = xf <- 'reference.tsv', sep = "\t", row.names = TRUE, col.names = NA)
  write.table(y, file = yf <- 'mixture.tsv', sep = "\t", row.names = TRUE, col.names = NA)
  message('OK')
  
  # run
  message("* Running CIBERSORT ... ", appendLF = FALSE)
  res <- CIBERSORT(xf, yf, ...)
  message('OK')
  
  stats <- setdiff(colnames(res), colnames(reference))
  list(coefficients = t(res[, !colnames(res) %in% stats]), stats = res[, stats])
}


#' Testing for Cell Type-Specific Gene Expression Differences
#' 
#' Uses `csSAM` (\cite{Shen-Orr2010}) to estimate cell type-specific 
#' expression and test for group differences in each cell type.
#' 
#' @param formula csSAM model to fit specified as a formula
#' @param data cell type proportion matrix (samples x cell types).
#' @param ... other arguments passed to [csSAM::csSAMfit].
#'  
#' @import csSAM
#' @export
bseqsc_csdiff <- function(formula, data = NULL, ...){
  
  csSAMfit(formula, data = data, ...)
  
}

#' Fitting EdgeR Model
#' 
#' @param x an [Biobase::ExpressionSet-class] object
#' @param formula EdgeR formula to fit
#' @param coef coefficient(s) of interest, for which a statistical test is performed.
#' 
#' @return an object of class `TopTags`, as returned by [edgeR::topTags].
#' 
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @export
fitEdgeR <- function(x, formula, coef){
  qlibrary('edgeR') # this is required for subsetting happening in glmQLFTest to work
  design <- model.matrix(formula, data = pData(x))
  esetA <- x[, rownames(design)]
  y <- DGEList(counts = exprs(esetA), genes = fData(esetA))
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = coef)
  top <- topTags(qlf, n = Inf)
  top$fit <- fit
  top
  
}
