%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotDensity.AffymetrixDataFile}
\alias{plotDensity.AffymetrixDataFile}
\alias{AffymetrixDataFile.plotDensity}
\alias{plotDensity.AffymetrixDataFile}
\alias{plotDensity,AffymetrixDataFile-method}


\title{Plots the densities of all samples}

\description{
 Plots the densities of all samples.
}

\usage{\method{plotDensity}{AffymetrixDataFile}(this, probes=1/20, types=NULL, ..., xlab=NULL, ylab="density (integrates to one)", log=TRUE, verbose=FALSE)}

\arguments{
  \item{probes}{The subset of probes to include.}
  \item{types}{The type of probes to include.}
  \item{...}{Additional arguments passed to
     \code{\link[aroma.affymetrix:getProbes.AffymetrixDataFile]{*getProbes}()}, \code{\link[aroma.affymetrix:getProbeIntensities.AffymetrixDataFile]{*getProbeIntensities}()} and
     \code{\link[aroma.light]{plotDensity.numeric}}.}
  \item{xlab,ylab}{The labels on the x and the y axes.}
  \item{log}{If \code{\link[base:logical]{TRUE}}, the density of the log (base 2) values are
     used, otherwise the non-logged values.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
