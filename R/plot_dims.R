#' Plot reduced dimensions
#'
#' Plot reduced dimensions from singlecellexperiment object and var them on color or plot expression level
#' @param sce_object A single cell experiment object
#' @param x Principal component to plot on x axis
#' @param y Principal component to plot on y axis
#' @param var Variable or gene name
#' @param palette Number of palette selected from qualitative series on colorbrewer2.org
#' @param hex_codes Manuel entry of HEX color codes to plot
#' @param point_size Size parameter of geom_point
#' @param alpha Decide the transparency of geom
#' @return Ggplot object plot of single-cell dataset with variables or genes as color scale
#' @export
plot_dims <- function(sce_object, x, y, var, palette = "Dark2", hex_codes = NA, point_size = 3, alpha = 1, theme = 18) {
    if (length(var) == sum(var %in% row.names(sce_object))) {
      if (length(var) > 1) {
        rowData <- t(data.frame(logcounts(sce_object)[var, ]))
        colData <- data.frame(Dim1 = sce_temp[[x]], Dim2 = sce_temp[[y]])
        temp <- cbind(rowData, colData)
        temp <- gather(temp, gene, logcounts, -Dim1, - Dim2)

        ggplot(temp, aes(Dim1, Dim2, col = logcounts)) +
          geom_point(size = point_size) +
          facet_wrap(~ gene) +
          #labs(x = x, y = y, color = group) +
          #guides(color = guide_colorbar(barwidth = 8, barheight = 1, ticks = FALSE, title.vjust = c(1.3), title = "Logcounts")) +
          scale_color_viridis(option = "inferno") +
          theme_bw(base_size = theme) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
      }
      else if (length(var) == 1) {
        temp <- data.frame(x.var = sce_object[[x]], y.var = sce_object[[y]], gene_name = Biobase::exprs(sce_object[var])[1, ])
        ggplot(temp, aes(x.var, y.var, col = gene_name), alpha = 1) +
          geom_point(size = point_size, alpha = alpha) +
          labs(x = x, y = y, color = var) +
          #guides(color = guide_colorbar(barwidth = 0.5, barheight = 8, ticks = FALSE)) +
          scale_color_viridis(option = "inferno") +
          theme_bw(base_size = theme) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      }
    }
    else if (length(var) != sum(var %in% row.names(sce_object))) {
      if (is.na(hex_codes) == TRUE) {
        temp <- data.frame(x.var = sce_object[[x]], y.var = sce_object[[y]], col = sce_object[[var]])
        ggplot(temp, aes(x.var, y.var, col = col), alpha = 0.8) +
          geom_point(size = point_size, alpha = alpha) +
          labs(x = x, y = y, color = var) +
          scale_color_brewer(type = "qual", palette = palette) +
          theme_bw(base_size = theme) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      }
      else {
        temp <- data.frame(x.var = sce_object[[x]], y.var = sce_object[[y]], col = sce_object[[var]])
        ggplot(temp, aes(x.var, y.var, col = col)) +
          geom_point(size = point_size, alpha = alpha) +
          labs(x = x, y = y, color = var) +
          scale_color_manual(values = hex_codes) +
          theme_bw(base_size = theme) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      }
    }
  }


