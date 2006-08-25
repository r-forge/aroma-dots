%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getProbesetIntensities.AffymetrixDataFile}
\alias{getProbesetIntensities.AffymetrixDataFile}
\alias{AffymetrixDataFile.getProbesetIntensities}
\alias{getProbesetIntensities.AffymetrixDataFile}
\alias{getProbesetIntensities,AffymetrixDataFile-method}


\title{Gets probe signals for a subset of probesets in one array}

\description{
  Gets probe signals for a subset of probesets in one array.
}

\usage{\method{getProbesetIntensities}{AffymetrixDataFile}(this, probesets, doTransform=isTransforming(this), ..., stratifyBy=c("nothing", "pmmm", "pm", "mm"), verbose=FALSE)}

\arguments{
 \item{probesets}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying probesets to be read.
   If \code{\link[base]{NULL}}, all probesets are read.}
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
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
