# Plotting functions
# 
# Author: Renaud Gaujoux
###############################################################################

#' @import ggplot2
.ggtheme <- function(x = NULL){
  
  # load color theme
  COLOR_THEME <- ldata('bseqscColorTheme')
  default <- theme_set(theme_bw(base_size = 12, base_family = "Helvetica"))
  theme_update(#panel.grid.major = element_line(size = 0.5, color = "grey"),
      axis.line = element_line(size = 0.7, color = "black")
      , legend.position = 'top'
      , legend.title = element_text(face = 'bold', size = 10)
      , text = element_text(size = 14)
      , legend.key = element_blank()
      , panel.grid.major.x = element_blank()
      , panel.grid.minor = element_blank()
      , strip.background = element_blank()
      , strip.text.x = element_text(face = 'bold')
      , strip.text.y = element_text(face = 'bold')
      , panel.border = element_rect(colour = "black"))
  
  COLOR_THEME <- sapply(unique(as.character(COLOR_THEME$item)), function(x) as.character(COLOR_THEME$color[COLOR_THEME$item %in% x]), simplify = FALSE)
  if( !is.null(x) ) invisible(COLOR_THEME[[x]])
  else invisible(COLOR_THEME)
}

#' Plots Total Count Per Cell Type
#' 
#' Generate a barplot
#' 
#' @param eset a matrix or `ExpressionSet` object.
#' @param cellType vector or name of phenotypic variable that contains the mapping 
#' of each column (i.e. each cell) to a cell type.
#' @param sampleID vector or name of phenotypic variable that contains the mapping 
#' of each column (i.e. each cell) to a biological sample.
#' 
#' @return a `ggplot` object with [ggplot2::geom_bar].
#' 
#' @export
#' @examples 
#' # generate random data
#' x <- rscData()
#' # check phenotypic variables
#' varLabels(x)
#' # plot counts
#' plotCellTotals(x, 'cellType', 'sampleID')
#' 
plotCellTotals <- function(eset, cellType, sampleID){
  
  df <- pData(eset)
  df$Total <- colSums(exprs(eset))
  mdf <- ddply(df, cellType, function(x){
        msample <- unlist(dlply(x, sampleID, function(x) mean(x$Total)))
        x$avgtot <- mean(msample)
        x$avgtot.sd <- sd(msample)
        x$lb <- x$avgtot - x$avgtot.sd
        x$ub <- x$avgtot + x$avgtot.sd
        x[1, ]
      })
  
  mdf[[cellType]] <- stats::reorder(mdf[[cellType]], -mdf$avgtot)
  
  # plot
  ggplot(mdf, aes_string(y = 'avgtot', x = cellType, fill = cellType)) + geom_bar(stat = 'identity') + 
      guides(fill = guide_legend('')) +  
      geom_errorbar(aes_string(ymin = 'lb', ymax = 'ub'), color = '#555555', width = .25) + 
      theme(axis.text.x = element_text(angle = 320, hjust = 0)) + ylab("Average total count") + xlab('')
  
}

#' Plotting Cell Type-Specific Expression Effect Size
#' 
#' Generates an horizontal barplot of cell type-specific expression 
#' differences.
#' 
#' @param x a `data.frame` as returned by [csSAM::csTopTable].
#' 
#' @return a `ggplot` object.
#' 
#' @importFrom dplyr mutate_
#' @export
plotEffectSize <- function(x){
  
  # reorder feature levels
  new_feat <- ~ stats::reorder(factor(Feature), rank(as.numeric(factor(Cell.type)) * max(abs(t)) + abs(t)))
  cst <- mutate_(x, .dots = setNames(list(new_feat), 'Feature'))
  
  # plot
  ggplot(cst, aes_string(x = 'Feature', y = 't', fill = 'Cell.type')) + 
      geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey30') + 
      ylab('Effect size') + xlab('Features') + 
      geom_bar(stat = 'identity', position = 'dodge') + coord_flip() 
  
  
}

#' P-Value Scatter Plot
#' 
#' @param formula Formula of the form `yvar ~ xvar` that specify which variable in `data`
#' contains variables for the x-axis and y-axis respectively.
#' @param data a `data.frame` object containing the data to plot
#' @param pval.th default p-value threshold used to filter te data shown in the plot.
#' Only data with \eqn{p <= pval.th} are used.
#' @param pval.th.y,pval.th.x p-value threshold to use for the variable in the x-axis and y-axis
#' respectively
#' @param label.th -log10(p-value) threshold above which data points are labelled.
#' Only the points in the lower triangle are labelled.
#' 
#' @export 
pvalueScatter <- function(formula, data, pval.th = 0.05, pval.th.y = pval.th, pval.th.x = pval.th, label.th = Inf){
  
  # parse formula
  t <- terms(formula, data = data)
  vars <- as.character(attr(t, 'variables'))[-1L]
  data$x <- data[[vars[[2L]]]]
  data$y <- data[[vars[[1L]]]]
  data <- data[data$x <= pval.th.x | data$y <= pval.th.y, , drop = FALSE]
  
  # add layer variables
  data$doLabel <- -log10(data$x) >= label.th & data$x <= data$y
  data$label <- ifelse(data$doLabel, data$Symbol, '')
  data$Improved <- data$Improved %||% data$x <= data$y
  
  # plot
  ggplot(data, aes_string(y = 'y', x = 'x', color = 'Improved')) + geom_point(size = .3) + geom_abline(slope = 1) + 
      scale_color_manual(values = c('TRUE' = '#aa2222', 'FALSE' = 'black'), guide = FALSE) + 
      scale_x_continuous(trans = reverselog_trans(10)) + scale_y_continuous(trans = reverselog_trans(10)) +
      geom_text(aes_string(label = 'label'), size = 2, hjust = -0.2, vjust = 1.2) + 
      xlab(vars[[2L]]) + ylab(vars[[1L]])

}

#' Generates Heatmap of Basis Matrix
#' 
#' Draws a heatmap of reference cell type expression profiles, highlighting 
#' marker genes.
#' 
#' @param x matrix or cell type reference expression profiles (genes x cell types)
#' @param markers named list of marker genes. The llist should be named using cell type names
#' corresponding to column names in `x`.
#' Genes are named using the identifiers corresponding to row names in `x`. 
#' Genes (rows) of `x` that are not in any element of `markers` are discarded.
#' Columns of `x` that match no element in `markers` are discarded.
#' @param scale scaling strategy. See [NMF::aheatmap].
#' @param ... other arguments passed to [NMF::aheatmap].
#' 
#' @importFrom NMF aheatmap
#' @importFrom AnnotationDbi unlist2
#' @export
plotBasis <- function(x, markers = NULL, scale = 'r1', ...){
  
  if( !is.null(markers) ){
    ml <- sapply(markers, intersect, rownames(x), simplify = FALSE)
    ml <- ml[lengths(ml) > 0]
    x <- x[unlist(ml), names(ml)]
    ann <- data.frame(Markers = factor(names(unlist2(ml)), levels = names(ml)))
  }
  aheatmap(x, scale = scale, annRow = ann, ...)
  
}

#' @importFrom scales trans_new log_breaks
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
      log_breaks(base = base), 
      domain = c(1e-100, Inf))
}
