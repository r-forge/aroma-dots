\name{plotBoxplot.AffymetrixCelSet}
\alias{plotBoxplot.AffymetrixCelSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to plot Boxplots }
\description{
    Create boxplots of values stored in AffymetrixCelSet object
}
\usage{
\method{plotBoxplot}{AffymetrixCelSet}(this, arrays = NULL, verbose = FALSE, calcArgs=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ an AffymetrixCelSet object }
  \item{arrays}{ arrays to calculate boxplot}
  \item{calcArgs}{named list of arguments to be passed to plotBoxplot.AffymetrixCelFile}
  \item{\dots}{ arguments passed to \code{bxp} for plotting boxplots }
}
\details{
The function calls plotBoxplot.AffymetrixCelFile for each array (without plotting the boxplot), so any changes in the calculation, such as the proportion of the IQR that the whiskers should go, etc., should be passed in this list.
}
\value{
    list appropriate to give as argument to \code{bxp}
}
\author{Elizabeth Purdom}
\seealso{ \code{\link{bxp}} }
