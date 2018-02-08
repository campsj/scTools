#' Select markers
#'
#' Select marker genes based on padj, auroc, cluster number and k and returns them as a vector
#' @param sce_object A single cell experiment object
#' @param k number of clusters
#' @param clust_num number of cluster to select
#' @param padj Adjusted p value
#' @param auroc Area Under Receiving Operating Curve
#' @param vector_name Name to give the returned vector
#' @return Vector of marker genes
#' @export

select_markers <- function(sce_object, k, clust_num, padj = 0.01, auroc = 0.85, vector_name = "genes") {
  #require(c("SC3", "lazyeval", "dplyr"))

  k_clusts <- paste0("sc3_", k, "_markers_clusts", sep = "")
  k_auroc <- paste0("sc3_", k, "_markers_auroc", sep = "")
  k_padj <- paste0("sc3_", k, "_markers_padj", sep = "")

  #sc3_plot_markers(sce_object, k = k, show_pdata = TRUE)

  genes <- sce_object %>%
    rowData() %>%
    as_tibble() %>%
    filter_(interp(~ a < x & b > y & c == z, a = as.name(k_padj), b = as.name(k_auroc), c = as.name(k_clusts), x = padj, y = auroc, z = clust_num)) %>%
    arrange_(k_padj)

  genes <- genes[ ,"mgi_symbol", drop = TRUE]
  assign(paste0(vector_name), genes, envir = .GlobalEnv)
}
