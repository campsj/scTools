#' Plot expression
#'
#' Plot variables or genes of interest on a variable number of principal components
#' @param sce_object A single cell experiment object
#' @param x Groups to plot on x axis
#' @param gene Principal component to plot on y axis
#' @param group Color variable
#' @param folder Name of designated folder to save image
#' @param violin_fill Fill color for violin
#' @param brew_type Kind color type to plot: sequential, diverging or qualitative
#' @param brew_palette Number of palette selected from qualitative series on colorbrewer2.org
#' @param point_size Size parameter of geom_point
#' @param scatter_width Measure how wide points should scatter
#' @param theme Number of theme_size
#' @param width Width of image in cm
#' @param height Height of image in cm
#' @return Gene expression plot
#' @export

plot_expression <- function(sce_object, x, gene, group, folder, violin_fill = "#f0f0f0", brew_type = "qual", brew_palette = "Set1", point_size = 1.5, scatter_width = 0.2, theme = 12,  width = 12, height = 10) {
  if (dir.exists(paste(folder, sep = "")) == FALSE)  {
    dir.create(paste(folder, sep =""))
  }
  gene %in% row.names(sce_object)
  temp <- data.frame(X = sce_object[[x]], Y = logcounts(sce_object[gene])[1, ], group = sce_object[[group]])
  ggplot(temp, aes(X, Y)) +
    geom_violin(fill = violin_fill) +
    geom_jitter(aes(col = group), width = scatter_width, size = point_size) +
    labs(x = x, y = paste0(gene, "Expression (logcounts)", sep = ""), color = group) +
    scale_color_brewer(type = brew_type, palette = brew_palette) +
    theme_bw(base_size = theme) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_text(size = 20, face = "bold")) +
    ggsave(paste0(folder, "/", gene, "_expression.tiff", sep = ""), width = width, height = height, units = "cm")
}
