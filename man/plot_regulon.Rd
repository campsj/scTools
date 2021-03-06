% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCENIC.R
\name{plot_regulon}
\alias{plot_regulon}
\title{Plot regulon activity on  PCA}
\usage{
plot_regulon(sce_object, regulon, binaryAct, group)
}
\arguments{
\item{sce_object}{SingleCellExperiment object}

\item{regulon}{Regulon name with genes in network between brackets}

\item{binaryAct}{Binary or AUC regulon information}

\item{group}{Grouping variable for shape}
}
\value{
Ggplot2 object
}
\description{
Plot regulon activity on  PCA
}
