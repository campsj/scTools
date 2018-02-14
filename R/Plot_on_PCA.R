#' Plot on PCA
#'
#' Plot variables or genes of interest on a variable number of principal components
#' @param sce_object A single cell experiment object
#' @param n Amount of principal component to plot
#' @param var Vector of variable names to plot
#' @return PCA plot of single-cell dataset with variables or genes as color scale
#' @export


plot_on_PCA <- function(sce_object, n, var) {
  for (v in var) {
  print(
    plotPCA(sce_object, ncomponents = n, colour_by = v) +
      ggtitle(v)
  )
  }
}
