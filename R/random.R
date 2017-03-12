# Generation of Random Single Cell Data
# 
# Author: Renaud Gaujoux
###############################################################################


#' Generates Random Single Cell Data
#' 
#' Generates random single cell data from a negative binomial distribution.
#' 
#' @param n number of features
#' @param p number of cells
#' @param celltypes number of cell types
#' @param samples number of biological samples
#' @param dispersion dispersion parameter passed to argument `size` of [stats::rnbinom].
#' @param .rng specification to set the random number with [rngtools::setRNG].
#' This enables to generate reproducible data.
#' 
#' @return an `ExpressionSet` object
#' 
#' @importFrom rngtools setRNG
#' @export
#' 
#' @examples 
#' x <- rscData()
#' x
#' varLabels(x)
#' summary(pData(x))
#' 
rscData <- function(n = 800, p = 5000, celltypes = 15, samples = 4, dispersion = .1, .rng = NULL){
  
  # setup RNG
  if( !is.null(.rng) ){
    orng <- setRNG(.rng)
    on.exit( setRNG(orng) )
  }
  
  # cell types
  ctnames <- paste0('CT_', seq(celltypes))
  ct <- gl(celltypes, p / celltypes, length = p, labels = ctnames)
  ct_table <- summary(ct)
  # samples
  sampleID <- sample(gl(samples, p / samples, length = p))
  # sample groups
  sgroup <- setNames(gl(2, samples / 2, length = samples, labels = c('A', 'B')), levels(sampleID))
  
  # generate from negative binomial distribution
  # sample mean of each cell type
  mu <- setNames(runif(celltypes, min = 1, max = 10), ctnames)
  cnames <- paste0('C', seq(p))
  x <- lapply(levels(ct), function(ctype){
        m <- mu[ctype]
        res <- rnbinom(n * ct_table[ctype], size = dispersion, mu = m)
        matrix(res, nrow = n)
      })
  x <- do.call(cbind, x)
  dimnames(x) <- list(paste0('G', seq(n)), cnames)
  pd <- data.frame(cellType = ct
                  , sampleID = sampleID # random assignment to biological samples
                  , group = sgroup[as.character(sampleID)]
                  , row.names = colnames(x))
  
  # wrap-up as an ExpressionSet object
  ExpressionSet(x, phenoData = AnnotatedDataFrame(pd))
  
}
