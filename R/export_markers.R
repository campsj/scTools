#' Filter and save markers
#'
#' Filters marker genes based on padj and k and saves them as a csv file
#' @param sce_object A single cell experiment object
#' @param k number of clusters
#' @param cluster number of cluster to select
#' @param padj Adjusted p value
#' @param auroc Area Under Receiving Operating Curve
#' @param vector_name Name to give the returned vector
#' @return Vector of marker genes
#' @export

export_markers <- function(sce_object, k, padj = 0.01, auroc = 0.8, folder, filename) {
  require(c("SC3", "lazyeval", "dplyr"))

  k_clusts <- paste0("sc3_", k, "_markers_clusts", sep = "")
  k_auroc <- paste0("sc3_", k, "_markers_auroc", sep = "")
  k_padj <- paste0("sc3_", k, "_markers_padj", sep = "")

  sc3_plot_markers(sce_object, k = k, show_pdata = TRUE)

  sce_object %>%
    rowData() %>%
    as_tibble() %>%
    group_by_(k_clusts) %>%
    filter_(interp(~ a < x & b > y, a = as.name(k_padj), b = as.name(k_auroc), x = padj, y = auroc)) %>%
    select_("mgi_symbol", "ensembl_gene_id", "log10_mean_counts", k_auroc, k_clusts, k_padj) %>%
    arrange_(k_padj, .by_group = TRUE) %>%
    write.csv(paste0(folder, "/", filename, ",.csv"))
}
