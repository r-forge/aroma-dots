\name{plotXvsY.AffymetrixCelFile}
\alias{plotXvsY.AffymetrixCelFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Scatterplots of AffymetrixCelFile objects }
\description{
    A function to plot scatterplot of two AffymetrixCelFile objects
}
\usage{
\method{plotXvsY}{AffymetrixCelFile}(this1, this2, indices = NULL, field="intensities", densityPlot = TRUE, loess = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this1}{ AffymetrixCelFile object (x axis)}
  \item{this2}{ AffymetrixCelFile object (y axis) }
  \item{indices}{ cell indices to include in plot; NULL will plot all }
  \item{field}{field to plot}
  \item{densityPlot}{ logical; whether to plot the points or just a density representation of them }
  \item{loess}{ logical; whether to plot loess curve}
  \item{\dots}{ arguments passed to plot or smoothScatter}
}
\details{
 Plots are done after taking log2 of the intensities field. Note that zero valued cells (in the original scale, not log scale) will not be plotted. 
}

\author{Elizabeth Purdom}
