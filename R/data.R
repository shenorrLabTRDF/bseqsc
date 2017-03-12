# Embedded Data
# 
# Author: Maayan Baron
# Author: Renaud Gaujoux
###############################################################################

#' Ribosomal and Mitochondrial Genes
#' 
#' @source http://www.genenames.org/cgi-bin/genefamilies/set/1054
#' 
#' @docType data
'MITRIB'

#' MITRIB Lists known ribosomal or mitochondrial genes.
#' 
#' @param identifier type of gene identifier:
#'   * 'SYMBOL': official gene symbols;
#'   * 'ENTREZID': NCBI Entrez gene IDs.
#' 
#' @return a character vector of gene IDs with names corresponding 
#' to their type (either 'RIB' or 'MIT').
#' 
#' @export
getMITRIB <- function(identifier = c('SYMBOL', 'ENTREZID')){
  
  identifier <- match.arg(identifier)
  # load mit-rib genes from package
  MITRIB <- ldata('MITRIB')
  setNames(as.character(MITRIB[[identifier]]), as.character(MITRIB$Type))
  
}


#' Marker Genes for Pancreas Islet Cells
#' 
#' List of marker genes for 6 pancreas cell sub-population:
#'   
#' @docType data
'pancreasMarkers'

#' Reference Basis Matrix of Pancreatic Islet Cell Types
#' 
#' Expression profiles for 6 pancreatic islet cell sub-populations, 
#' over a set of marker genes that are characteristic of each cell type.
#' It is used to estimate cell type proportion from bulk tissue total gene 
#' expression.
#' 
#' @docType data
'PancreasIslet'

#' Color Theme for BSEQsc Vignettes
#' 
#' @docType data
'bseqscColorTheme'
