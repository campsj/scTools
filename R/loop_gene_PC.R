#' Loop over genes and components vector
#'
#' Plot variables or genes of interest on a variable number of principal components
#' @param genes Vector of genes names
#' @param pcs List of principals components to plot
#' @param dir Name of directory to save files
#' @param subdir Name of subdirectory to save files
#' @param sce_object A single cell experiment object
#' @param alpha Decide the transparency of geom
#' @param width Width of image in cm
#' @param height Height of image in cm
#' @return PCA plot of single-cell dataset with expression of selected genes as color scale
#' @export

loop_gene_PC <- function(genes, pcs, dir = dir, subdir = subdir, sce_object = sce_object, alpha = alpha, width = "14", height = "10", units = "cm") {
if (dir.exists(paste(dir, "/", subdir, sep = "")) == FALSE)  {
  dir.create(paste(dir, "/", subdir, sep =""))
}
if (exists("pcs") == FALSE) {
  print("Create pcs: list with principal components")
}
for(g in genes) {
  if(g %in% row.names(sce_object)){
    for (i in 1:6) {
      plot_components(sce_object, pcs[i, 1], pcs[i, 2], group = g, folder = subdir, alpha = 0.8, gene = TRUE, width = width, height = height, units = units)
    }
  }
}
}
