#' Plot regulon activity on  PCA
#'
#' @param sce_object SingleCellExperiment object
#' @param regulon Regulon name with genes in network between brackets
#' @param binaryAct Binary or AUC regulon information
#' @param group Grouping variable for shape
#' @return Ggplot2 object
#' @export

plot_regulon <- function(sce_object, regulon, binaryAct, group) {
  label <- unlist(strsplit(regulon, " "))[1]
  label <- unlist(strsplit(label, "_"))[1]
  df <- data.frame(PC1 = sce_object[["PC1"]], PC2 = sce_object[["PC2"]], regulon_activity = binaryAct[regulon, colnames(sce_object)],
                   group = sce_object[[group]])
  ggplot2::ggplot(df, aes(x = PC1, y = PC2, col = regulon_activity)) +
    geom_point(size = 3) +
    labs(title = label, x = "Principal component 1", y = "Principal component 2", color = "Regulon\nactivity") +
    guides(shape = guide_legend(title = element_blank())) +
    viridis::scale_color_viridis(option = "plasma", guide = guide_colorbar(ticks = FALSE, barwidth = 0.5,
                                                                           barheight = 3, title.vjust = 0.9,
                                                                           title.position = "top")) +
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "right", legend.title = element_text(angle = 0), axis.title = element_text(size = 16),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5), legend.justification = c(0, 1))
}
