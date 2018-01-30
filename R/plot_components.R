#' Plot components
#'
#' Plot variables or genes of interest on a variable number of principal components
#' @param sce_object A single cell experiment object
#' @param PCx Principal component to plot on x axis
#' @param PCy Principal component to plot on y axis
#' @param group Variable or gene name
#' @param folder Name of designated folder to save images
#' @param alpha Decide the transparency of geom
#' @param palette Number of palette selected from qualitative series on colorbrewer2.org
#' @param gene If TRUE it will plot genes and no variables
#' @param width Width of image in cm
#' @param height Height of image in cm
#' @return PCA plot of single-cell dataset with variables or genes as color scale
#' @export
plot_components <- function(sce_object, PCx, PCy, group, folder, alpha = 0.8, palette = 2, gene = FALSE, width = 14, height = 10, units = "cm") {
  if (group == "genotype") {
    temp <- data.frame(PCa = sce_object[[PCx]], PCb = sce_object[[PCy]], col = sce_object[[group]])
    ggplot(temp, aes(PCa, PCb, col = col)) +
      geom_point(size = 3, alpha = alpha) +
      labs(x = paste("Principal component ", PCx, sep = ""), y = paste("Principal compoenent ", PCy, sep = ""), color = paste(group)) +
      #guides() +
      scale_color_manual(values = c("#FFC90C", "#89CAFF")) +
      theme_bw(base_size=18) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_blank()) +
      ggsave(paste("plots/",folder,"/", group,"_", PCx, "_", PCy,".tiff", sep=""), width = width, height = height, units = units)
  } else if (group == g) {
    temp <- data.frame(PCa = sce_object[[PCx]], PCb = sce_object[[PCy]], gene_name = exprs(sce_object[group])[1, ])
    ggplot(temp, aes(PCa, PCb, col = gene_name), alpha = 0.8) +
      geom_point(size = 3, alpha = alpha) +
      labs(x = paste("Principal component ", PCx, sep = ""), y = paste("Principal compoenent ", PCy, sep = ""), color = paste(group)) +
      guides(color = guide_colorbar(barwidth = 0.5, barheight = 8, ticks = FALSE)) +
      scale_color_viridis(option = "viridis") +
      theme_bw(base_size=18) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_text(size = 24, face = "bold")) +
      ggsave(paste("plots/",folder,"/",group,"_", PCx, "_", PCy,".tiff", sep=""), width = width, height = height, units = units)
  } else {
    temp <- data.frame(PCa = sce_object[[PCx]], PCb = sce_object[[PCy]], col = sce_object[[group]])
    ggplot(temp, aes(PCa, PCb, col = col), alpha = 0.8) +
      geom_point(size = 3, alpha = alpha) +
      labs(x = paste("Principal component ", PCx, sep = ""), y = paste("Principal compoenent ", PCy, sep = ""), color = paste(group)) +
      #guides() +
      scale_color_brewer(type = "qual", palette = palette) +
      theme_bw(base_size=18) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_blank()) +
      ggsave(paste("sc3/",folder,"/",group,"_", PCx, "_", PCy,".tiff", sep=""), width = width, height = height, units = units)
  }
}
