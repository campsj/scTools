#' Plot expression
#'
#' Plot variables or genes of interest on a variable number of principal components
#' @param sce_object A single cell experiment object
#' @param var Gene or genes to plot
#' @param group Grouping variable
#' @param theme Size of theme
#' @return Gene expression plot
#' @export

plot_expression <- function(sce_object, var, group, theme = 12) {

  rowData <- t(data.frame(logcounts = logcounts(sce_object)[var, ]))
  colData <- data.frame(cluster = sce_object[[group]])
  temp <- cbind(rowData, colData)
  temp <- gather(temp, gene, logcounts, -cluster)

  ggplot(temp, aes(cluster, logcounts)) +
    geom_violin(aes(fill = cluster), scale = "width") +
    #geom_jitter(aes(col = genotype), width = 0.2, size = 2) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    facet_grid(. ~ gene, scales = "free") +
    labs(y = "Expression (logcounts)", fill = "Cluster") +
    #scale_fill_manual(values = c("#1f78b4", "#b2df8a")) +
    scale_fill_brewer() +
    theme_bw(base_size = theme) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_text(size = 16, face = "bold"), legend.position = "top", axis.title.x = element_blank(),
          strip.text.x = element_text(size = 14, face = "bold.italic", angle = 0), axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  }


