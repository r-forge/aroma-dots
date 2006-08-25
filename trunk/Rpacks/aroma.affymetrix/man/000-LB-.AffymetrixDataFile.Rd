%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[.AffymetrixDataFile}
\alias{[.AffymetrixDataFile}
\alias{AffymetrixDataFile.[}
\alias{[.AffymetrixDataFile}
\alias{[,AffymetrixDataFile-method}


\title{Gets the probe intensities for a set of probes}

\description{
  Gets the probe intensities for a set of probes.
}

\usage{\method{[}{AffymetrixDataFile}(this, i=NULL, drop=FALSE)}

\arguments{
 \item{i}{An index \code{\link[base]{vector}} specifying which probe to be queried.
   If \code{\link[base]{NULL}}, all probes are considered.}
 \item{drop}{If \code{\link[base:logical]{TRUE}} and the length of the return vector is less or
   equal to one, then its dimensions are dropped, otherwise not.}
}

\value{
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getProbeIntensities.AffymetrixDataFile]{*getProbeIntensities}()}.
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
