#' Plot reduced dimensions
#'
#' Plot reduced dimensions from singlecellexperiment object and var them on color or plot expression level
#' @param sce_object A single cell experiment object
#' @param x Principal component to plot on x axis
#' @param y Principal component to plot on y axis
#' @param color factor for color
#' @param shape factor for shape
#' @param labels list containing labels for gene names
#' @param col_values Manuel entry of colors or HEX color codes to plot
#' @param point_size Size parameter of geom_point
#' @param alpha Decide the transparency of geom
#' @param theme numerical specifying base size of theme
#' @param nrow numerical specifying the amount of rows
#' @param ncol numerical specifying the amount of columns
#' @return Ggplot object plot of single-cell dataset with variables or genes as color scale
#' @export
plot_dims <- function(sce_object, x = "PC1", y = "PC2", color, shape = NA, labels = NA, col_values = NA, point_size = 3, alpha = 1, theme = 12, nrow = NULL, ncol = NULL) {
    if (sum(color %in% row.names(sce_object)) >= 1) {
      if (length(color) > 1) {
        rowData <- NULL
        #rowData <- t(data.frame(logcounts(sce_object)[color, ]))
        for (gene in color) {
          if (gene %in% row.names(sce_object)) {
            rowData[[gene]]  <- logcounts(sce_object)[gene, ]
          }
          else {
            print(paste0(gene, " not expressed or written wrong!", sep = ""))
          }
        }
        rowData <- data.frame(rowData)
        colData <- data.frame(Dim1 = sce_object[[x]], Dim2 = sce_object[[y]])
        temp <- cbind(rowData, colData)
        temp <- tidyr::gather(temp, gene, logcounts, -Dim1, - Dim2)

        ggplot(temp, aes(Dim1, Dim2, col = logcounts)) +
          geom_point(size = point_size) +
          facet_wrap(~ gene, labeller = labeller(gene = labels), nrow = nrow, ncol = ncol) +
          #labs(x = x, y = y, color = group) +
          #guides(color = guide_colorbar(barwidth = 8, barheight = 1, ticks = FALSE, title.vjust = c(1.3), title = "Logcounts")) +
          viridis::scale_color_viridis(option = "inferno", guide = guide_colourbar(ticks = FALSE)) +
          theme_bw(base_size = theme) +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(),
                axis.title = element_blank(), axis.line = element_blank(), strip.background = element_blank(),
                strip.text.x = element_text(face = "bold.italic", size = theme + 4), legend.justification = c(0, 1))
      }
      else if (length(color) == 1) {
        temp <- data.frame(x.var = sce_object[[x]], y.var = sce_object[[y]], gene_name = Biobase::exprs(sce_object[color])[1, ])
        ggplot(temp, aes(x.var, y.var, col = gene_name), alpha = 1) +
          geom_point(size = point_size, alpha = alpha) +
          labs(x = x, y = y, color = color) +
          #guides(color = guide_colorbar(barwidth = 0.5, barheight = 8, ticks = FALSE)) +
          viridis::scale_color_viridis(option = "inferno") +
          theme_bw(base_size = theme) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      }
    }
    else if (sum(color %in% row.names(sce_object)) == 0) {
      if (is.na(shape) == FALSE) {
        temp <- data.frame(x.var = sce_object[[x]], y.var = sce_object[[y]], col.var = sce_object[[color]], shape.var = sce_object[[shape]])
      }
      else {
        shape.var <- c("1")
        temp <- data.frame(x.var = sce_object[[x]], y.var = sce_object[[y]], col.var = sce_object[[color]])
      }
      if (is.na(col_values) == TRUE) {
        ggplot(temp, aes(x.var, y.var, col = col.var, shape = shape.var)) +
          geom_point(size = point_size, alpha = alpha) +
          labs(x = x, y = y, color = color) +
          theme_bw(base_size = theme) +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), axis.text = element_blank(), axis.ticks = element_blank())
      }
      else {
        ggplot(temp, aes(x.var, y.var, col = col.var, shape = shape.var)) +
          geom_point(size = point_size, alpha = alpha) +
          labs(x = x, y = y, color = color) +
          scale_color_manual(values = col_values) +
          theme_bw(base_size = theme) +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), axis.text = element_blank(), axis.ticks = element_blank())
      }
    }
  }


