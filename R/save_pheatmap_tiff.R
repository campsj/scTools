#' Save heatmap
#'
#' Save heatmap in designated folder
#' @param x Object
#' @param filename Name of output file
#' @param width Width of image in cm
#' @param height Height of image in cm
#' @return PCA plot of single-cell dataset with variables or genes as color scale
#' @export

save_pheatmap_tiff <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  tiff(filename, width=width, height=height, units = "cm", res = 300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
