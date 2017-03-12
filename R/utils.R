# Utility functions
# 
# Author: Renaud Gaujoux
# Copyright Renaud Gaujoux (2016)
###############################################################################

#' @import xbioc pkgmaker stringr plyr
#' @bibliography ~/Documents/articles/library.bib
"_PACKAGE"

.PACKAGE_NAME <- 'bseqsc'

`%||%` <- function(a, b) if( is.null(a) ) b else a

#' Setup BSeq-SC External Dependencies
#' 
#' Configure `bseqsc` by setting up `CIBERSORT` source code.
#' 
#' `BSeq-sc` uses `CIBERSORT` to estimate cell type proportions, based on reference expression profiles.
#'  Due to licensing requirements, source code for this algorithm needs to be
#'  downloaded separately from its website http://cibersort.stanford.edu.
#'  It is released under the Stanford Non-Commercial License.
#' In order to use it with `bseqsc` you will need to:
#' 
#'   1. Got to http://cibersort.stanford.edu
#'   2. Register and log in
#'   3. Download the latest R source code from the [Download
#'      section](http://cibersort.stanford.edu/download.php). 
#'   4. Configure `bseqsc` by pointing it to the downloaded file. This is
#'      done using the function `bseqsc_config`, which will copy the given R
#'      source file into the `R-data/bseqsc` sub-directory in the user's home
#'      directory for subsequent usage:
#' 
#' ```
#' bseqsc_config('path/to/downloaded/source/CIBERSORT.R')
#' ```
#' 
#' @param file path to the CIBERSORT source R code.
#' @param error logical that indicates if an error should be thrown 
#' if configuration failed.
#' 
#' @return the path where the file was copied, or `NULL` if `bseqsc` is not correctly 
#' configured.
#' 
#' @export
bseqsc_config <- function(file = NULL, error = FALSE){
	
	
	path <- userData(package = .PACKAGE_NAME, create = TRUE)
	cib <- file.path(path, 'CIBERSORT.R')
	if( !is.null(file) ){
		f <- file
		if( length(f) )
			file.copy(f, cib, overwrite = TRUE)
	}
	
	if( !file.exists(cib) ){
		if( !error ) return(NULL)
		stop("Could not find CIBERSORT source code at '", dirname(cib), "'.  Please ensure you correctly configured bseqsc. See ?bseqsc_config.\n")
	}
	
	env <- parent.frame()
	source(cib, local = env)
  # return main function
  res <- getFunction('CIBERSORT', where = env)
  attr(res, 'file') <- cib
	invisible(res)
}

count2tpm <- function(x, C = 10^6){
  tot <- colSums(x)
  sweep(x, 2L, tot, '/') * C
  
}
