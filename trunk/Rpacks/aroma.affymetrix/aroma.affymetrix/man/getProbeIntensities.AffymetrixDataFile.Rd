%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getProbeIntensities.AffymetrixDataFile}
\alias{getProbeIntensities.AffymetrixDataFile}
\alias{AffymetrixDataFile.getProbeIntensities}
\alias{getProbeIntensities.AffymetrixDataFile}
\alias{getProbeIntensities,AffymetrixDataFile-method}


\title{Gets the probe intensities for a set of probes}

\description{
  Gets the probe intensities for a set of probes.
}

\usage{\method{getProbeIntensities}{AffymetrixDataFile}(this, probes=NULL, ..., verbose=FALSE, .readMap=getReadMap(this), .doTransform=isTransforming(this))}

\arguments{
 \item{probes}{Am index \code{\link[base]{vector}} specifying which probe to be queried.
   If \code{\link[base]{NULL}}, all probes are considered.}
 \item{...}{Not used.}
 \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
 \item{.readMap}{(Internal) A read map}
 \item{.doTransform}{(Internal) Transform signals or not?}
}

\value{
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
