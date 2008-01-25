\name{plotBoxplot.AffymetrixCelFile}
\alias{plotBoxplot.AffymetrixCelFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to plot boxplots }
\description{
  Create boxplot of values stored in AffymetrixCelFile object
}
\usage{
\method{plotBoxplot}{AffymetrixCelFile}(this, subset = 1/2, verbose = FALSE, field = "intensities", trans = log2, types = "pm", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ an AffymetrixCelFile object }
  \item{subset}{ passed to \code{identifyCells} to use only a subset of probes to calculate boxplot }
  \item{field}{ name of field to create boxplot of }
  \item{trans}{ function describing a transformation wanted of the data (must take vectors)}
  \item{types}{ passed to \code{identifyCells} }
  \item{\dots}{ arguments passed to \code{bxp} for plotting boxplots }
}
\value{
    list appropriate to give as argument to \code{bxp}
}
\author{Elizabeth Purdom}
\seealso{ \code{\link{bxp}} }
