%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getProbeIntensities.AffymetrixDataSet}
\alias{getProbeIntensities.AffymetrixDataSet}
\alias{AffymetrixDataSet.getProbeIntensities}
\alias{getProbeIntensities.AffymetrixDataSet}
\alias{getProbeIntensities,AffymetrixDataSet-method}


\title{Gets probe intensities from a set of probes and a set of arrays}

\description{
  Gets probe intensities from a set of probes and a set of arrays.
}

\usage{\method{getProbeIntensities}{AffymetrixDataSet}(this, arrays=NULL, probes=NULL, doTransform=isTransforming(this), verbose=FALSE, ...)}

\arguments{
 \item{arrays}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying arrays to be read.
   If \code{\link[base]{NULL}}, all arrays are read.}
 \item{probes}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying probes to be read.
   If \code{\link[base]{NULL}}, all probes are read.}
 \item{doTransform}{If \code{\link[base:logical]{TRUE}}, probe signals are transformed by the
   internal transformation functions, otherwise not.}
 \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
   May also be a \code{\link[R.utils]{Verbose}} object.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{numeric}} \eqn{NxK} matrix, where \eqn{N} is the number of
  probes read, and \eqn{K} is the number of arrays read.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
