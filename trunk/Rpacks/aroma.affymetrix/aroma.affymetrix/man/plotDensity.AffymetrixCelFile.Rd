%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotDensity.AffymetrixCelFile}
\alias{plotDensity.AffymetrixCelFile}
\alias{AffymetrixCelFile.plotDensity}
\alias{plotDensity.AffymetrixCelFile}
\alias{plotDensity,AffymetrixCelFile-method}


\title{Plots the density of the probe signals on the array}

\description{
 Plots the density of the probe signals on the array.
}

\usage{\method{plotDensity}{AffymetrixCelFile}(this, subset=1/2, types=NULL, ..., xlab=NULL, ylab="density (integrates to one)", log=TRUE, annotate=TRUE, verbose=FALSE)}

\arguments{
  \item{subset}{The subset of probes to include.}
  \item{types}{The type of probes to include.}
  \item{...}{Additional arguments passed to
     \code{identifyCells()} of \code{\link{AffymetrixCdfFile}},
     , \code{\link[aroma.affymetrix:getData.AffymetrixCelFile]{*getData}()}, and
     \code{\link[aroma.light]{plotDensity.numeric}}.}
  \item{xlab,ylab}{The labels on the x and the y axes.}
  \item{log}{If \code{\link[base:logical]{TRUE}}, the density of the log (base 2) values are
     used, otherwise the non-logged values.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson \url{http://www.braju.com/R/}}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
