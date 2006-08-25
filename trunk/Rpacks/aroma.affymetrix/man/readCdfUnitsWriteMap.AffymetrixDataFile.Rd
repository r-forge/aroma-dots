%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.MAP.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readCdfUnitsWriteMap.AffymetrixDataFile}
\alias{readCdfUnitsWriteMap.AffymetrixDataFile}
\alias{AffymetrixDataFile.readCdfUnitsWriteMap}
\alias{readCdfUnitsWriteMap.AffymetrixDataFile}
\alias{readCdfUnitsWriteMap,AffymetrixDataFile-method}


\title{Read the probe map for a subset of units}

\description{
  Read the probe map for a subset of units.
}

\usage{\method{readCdfUnitsWriteMap}{AffymetrixDataFile}(this, units=NULL, ...)}

\arguments{
 \item{units}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of indices of units to be read.
   If \code{\link[base]{NULL}}, all units are considered.}
 \item{...}{Not used.}
}

\value{
  Returns an \code{\link[base]{integer}} \code{\link[base]{vector}} which is a \emph{write} map.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[affxparser]{readCdfUnitsWriteMap}} is used.
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
