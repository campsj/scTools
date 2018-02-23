#' Plot expression of marker genes identified by \code{SC3} as a heatmap.
#'
#' By default the genes with the area under the ROC curve (AUROC) > 0.85
#' and with the p-value < 0.01 are selected and the top 10 marker
#' genes of each cluster are visualized in this heatmap.
#'
#' @name plot_markers
#' @aliases sc3_plot_markers, sc3_plot_markers,SingleCellExperiment-method
#'
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param auroc area under the ROC curve
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#'
#' @importFrom pheatmap pheatmap
plot_markers.SingleCellExperiment <- function(object, k, auroc, p.val, show_pdata) {
  if (is.null(metadata(object)$sc3$consensus)) {
    warning(paste0("Please run sc3_consensus() first!"))
    return(object)
  }
  hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
  dataset <- get_processed_dataset(object)
  if (!is.null(metadata(object)$sc3$svm_train_inds)) {
    dataset <- dataset[, metadata(object)$sc3$svm_train_inds]
  }

  add_ann_col <- FALSE
  ann <- NULL
  if (!is.null(show_pdata)) {
    ann <- make_col_ann_for_heatmaps(object, show_pdata)
    if (!is.null(ann)) {
      add_ann_col <- TRUE
      # make same names for the annotation table
      rownames(ann) <- colnames(dataset)
    }
  }

  # get all marker genes
  markers <- organise_marker_genes(object, k, p.val, auroc)
  # get top 10 marker genes of each cluster
  markers <- markers_for_heatmap(markers)

  row.ann <- data.frame(Cluster = factor(markers[, 1], levels = unique(markers[, 1])))
  rownames(row.ann) <- markers$feature_symbol

  do.call(pheatmap::pheatmap, c(list(dataset[markers$feature_symbol, , drop = FALSE], show_colnames = FALSE,
                                     cluster_rows = FALSE, cluster_cols = hc, cutree_cols = k, annotation_row = row.ann, annotation_names_row = FALSE,
                                     gaps_row = which(diff(markers[, 1]) != 0), cellheight = 10), list(annotation_col = ann)[add_ann_col]))
}

#' @rdname plot_markers
#' @aliases plot_markers
setMethod("sc3_plot_markers", signature(object = "SingleCellExperiment"), sc3_plot_markers.SingleCellExperiment)
