#' Plot components
#'
#' Plot variables or genes of interest on a variable number of principal components
#' @param sce_object A single cell experiment object
#' @param PCx Principal component to plot on x axis
#' @param PCy Principal component to plot on y axis
#' @param group Variable or gene name
#' @param folder Name of designated folder to save image
#' @param subfolder Name of designated subfolder to save image
#' @param gene If TRUE it will plot genes and no variables
#' @param brewer If TRUE it will use colorpalettes from colorbrewer2.org
#' @param palette Number of palette selected from qualitative series on colorbrewer2.org
#' @param hex_codes Manuel entry of HEX color codes to plot
#' @param point_size Size paramter of geom_point
#' @param alpha Decide the transparency of geom
#' @param width Width of image in cm
#' @param height Height of image in cm
#' @param units Units to plot image in
#' @return PCA plot of single-cell dataset with variables or genes as color scale
#' @export
plot_components <- function(sce_object, PCx, PCy, group, folder, subfolder, gene = FALSE, brewer = TRUE, palette = 2, hex_codes, point_size = 3, alpha = 0.8, theme = 18, width = 14, height = 10, units = "cm") {
  if (dir.exists(paste(folder, "/", subfolder, sep = "")) == FALSE)  {
    dir.create(paste(folder, "/", subfolder, sep =""))
  }
  if (brewer == FALSE) {
    temp <- data.frame(PCa = sce_object[[PCx]], PCb = sce_object[[PCy]], col = sce_object[[group]])
    ggplot(temp, aes(PCa, PCb, col = col)) +
      geom_point(size = point_size, alpha = alpha) +
      labs(x = paste("Principal component ", PCx, sep = ""), y = paste("Principal compoenent ", PCy, sep = ""), color = paste(group)) +
      scale_color_manual(values = hex_codes) +
      theme_bw(base_size = theme) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_blank()) +
      ggsave(paste(folder, "/",subfolder, "/", group, "_", PCx, "_", PCy,".tiff", sep=""), width = width, height = height, units = units)
  } else if (gene == TRUE) {
    temp <- data.frame(PCa = sce_object[[PCx]], PCb = sce_object[[PCy]], gene_name = exprs(sce_object[group])[1, ])
    ggplot(temp, aes(PCa, PCb, col = gene_name), alpha = 0.8) +
      geom_point(size = point_size, alpha = alpha) +
      labs(x = paste("Principal component ", PCx, sep = ""), y = paste("Principal compoenent ", PCy, sep = ""), color = paste(group)) +
      guides(color = guide_colorbar(barwidth = 0.5, barheight = 8, ticks = FALSE)) +
      scale_color_viridis(option = "viridis") +
      theme_bw(base_size = theme) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_text(size = 24, face = "bold")) +
      ggsave(paste(folder, "/", subfolder, "/", group, "_", PCx, "_", PCy,".tiff", sep=""), width = width, height = height, units = units)
  } else {
    temp <- data.frame(PCa = sce_object[[PCx]], PCb = sce_object[[PCy]], col = sce_object[[group]])
    ggplot(temp, aes(PCa, PCb, col = col), alpha = 0.8) +
      geom_point(size = point_size, alpha = alpha) +
      labs(x = paste("Principal component ", PCx, sep = ""), y = paste("Principal compoenent ", PCy, sep = ""), color = paste(group)) +
      scale_color_brewer(type = "qual", palette = palette) +
      theme_bw(base_size = theme) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_blank()) +
      ggsave(paste(folder, "/", subfolder, "/", group, "_", PCx, "_", PCy,".tiff", sep=""), width = width, height = height, units = units)
  }
}
