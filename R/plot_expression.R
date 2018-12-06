#' Plot expression
#'
#' Plot variables or genes of interest on a variable number of principal components
#' @param sce_object A single cell experiment object
#' @param var Gene or genes to plot
#' @param group Grouping variable
#' @param theme Size of theme
#' @param facet Orientation of facets
#' @param type type of plot, option for violin or boxplot
#' @param scale Set scale of violin, if "area" (default), all violins have the same area (before trimming the tails). If "count", areas are scaled proportionally to the number of observations. If "width", all violins have the same maximum width.
#' @return Gene expression plot
#' @import tidyr
#' @export

plot_expression <- function(sce_object, var, group, theme = 12, facet = "vertical", type = "violin", scale = "width") {
  rowData <- NULL
  for (gene in var) {
    if (gene %in% row.names(sce_object)) {
      rowData[[gene]]  <- logcounts(sce_object)[gene, ]
    }
    else {
      print(paste0(gene, " not expressed or written wrong!", sep = ""))
    }

  rowData <- data.frame(rowData)
  colData <- data.frame(cluster = sce_object[[group]])
  temp <- cbind(rowData, colData)
  temp <- tidyr::gather(temp, gene, logcounts, -cluster)
  temp$gene <- factor(temp$gene, levels = var)

  p <- ggplot(temp, aes(cluster, logcounts))
  }
  if (type == "violin") {
    p <- p + geom_violin(aes(fill = cluster, col = cluster), scale = scale)
  }
  else if (type == "boxplot") {
    p <- p + geom_boxplot(aes(fill = cluster))
  }


if (facet == "vertical") {
  p +
    facet_grid(gene ~ ., scales = "free") +
    labs(y = "Expression (logcounts)", fill = "Cluster") +
    #scale_fill_manual(values = c("#1f78b4", "#b2df8a")) +
    #scale_fill_brewer() +
    theme_bw(base_size = theme)
}
else if (facet == "horizontal") {
  p +
    facet_grid(. ~ gene, scales = "free") +
    labs(y = "Expression (logcounts)", fill = "Cluster") +
    #scale_fill_manual(values = c("#1f78b4", "#b2df8a")) +
    #scale_fill_brewer() +
    theme_bw(base_size = theme)
}
  }


