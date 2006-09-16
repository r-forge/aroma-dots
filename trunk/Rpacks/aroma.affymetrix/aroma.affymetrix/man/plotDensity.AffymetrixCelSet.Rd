%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotDensity.AffymetrixCelSet}
\alias{plotDensity.AffymetrixCelSet}
\alias{AffymetrixCelSet.plotDensity}
\alias{plotDensity.AffymetrixCelSet}
\alias{plotDensity,AffymetrixCelSet-method}


\title{Plots the densities of all samples}

\description{
 Plots the densities of all samples.
}

\usage{\method{plotDensity}{AffymetrixCelSet}(this, subset=1/2, types=NULL, ..., col=seq(this), lty=NULL, lwd=NULL, annotate=TRUE, add=FALSE, verbose=FALSE)}

\arguments{
  \item{subset}{The subset of probes to include.}
  \item{types}{The type of probes to include.}
  \item{...}{Additional arguments passed to
     \code{\link{plotDensity.AffymetrixDataFile}}.}
  \item{col}{A \code{\link[base]{vector}} of colors for each of the arrays.}
  \item{lty}{A \code{\link[base]{vector}} of line types for each of the arrays.}
  \item{lwd}{A \code{\link[base]{vector}} of line widths for each of the arrays.}
  \item{add}{If \code{\link[base:logical]{FALSE}}, a new plot is created, otherwise the generated
    graphics is added to the current plot.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
