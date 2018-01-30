#' Extract principal component gene weights
#'
#' Save heatmap in designated folder
#' @param sce_object A single cell experiment object
#' @param ntop Amount of genes to use for PCA
#' @param file_name Name for output csv file
#' @return Extract PCA gene weights in dataframe pca_loadings and save as csv
#' @export

extract_PC_gene_weights <- function(sce_object, ntop = 500, file_name = "PCA_gene_weights.csv") {
  rv <- rowVars(assay(sce_object))
  selected_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(sce_object)[selected_genes,]))
  pca_loadings <- as.data.frame(pca$rotation)
  pca_loadings[order(pca_loadings[,1], decreasing = TRUE),]
  write.csv(pca_loadings, file_name)
}
