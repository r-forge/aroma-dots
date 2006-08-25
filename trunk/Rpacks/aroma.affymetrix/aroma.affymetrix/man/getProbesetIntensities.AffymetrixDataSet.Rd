%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getProbesetIntensities.AffymetrixDataSet}
\alias{getProbesetIntensities.AffymetrixDataSet}
\alias{AffymetrixDataSet.getProbesetIntensities}
\alias{getProbesetIntensities.AffymetrixDataSet}
\alias{getProbesetIntensities,AffymetrixDataSet-method}


\title{Gets probe signals for a subset of probesets and a subset of arrays}

\description{
  Gets probe signals for a subset of probesets and a subset of arrays.
}

\usage{\method{getProbesetIntensities}{AffymetrixDataSet}(this, probesets=NULL, readMap=getReadMap(this), doTransform=isTransforming(this), ..., verbose=FALSE)}

\arguments{
 \item{arrays}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying arrays to be read.
   If \code{\link[base]{NULL}}, all arrays are read.}
 \item{probesets}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying probesets to be read.
   If \code{\link[base]{NULL}}, all probesets are read.}
 \item{readMap}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying a read map used to read the
   probe signals.}
 \item{doTransform}{If \code{\link[base:logical]{TRUE}}, probe signals are transformed by the
   internal transformation functions, otherwise not.}
 \item{...}{Arguments passed to the low-level function for read probesets,
   e.g. \code{\link[affxparser]{readCelUnits}} or \code{\link[aroma.apd]{readApdUnits}}.}
 \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
   May also be a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns a named \code{\link[base]{list}} structure.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
