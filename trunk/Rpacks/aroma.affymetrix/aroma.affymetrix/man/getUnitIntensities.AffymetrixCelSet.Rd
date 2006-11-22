%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitIntensities.AffymetrixCelSet}
\alias{getUnitIntensities.AffymetrixCelSet}
\alias{AffymetrixCelSet.getUnitIntensities}
\alias{getUnitIntensities.AffymetrixCelSet}
\alias{getUnitIntensities,AffymetrixCelSet-method}


\title{Gets cell signals for a subset of units and a subset of arrays}

\description{
  Gets cell signals for a subset of units and a subset of arrays.
}

\usage{\method{getUnitIntensities}{AffymetrixCelSet}(this, units=NULL, ..., verbose=FALSE)}

\arguments{
 \item{arrays}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying arrays to be read.
   If \code{\link[base]{NULL}}, all arrays are read.}
 \item{units}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying units to be read.
   If \code{\link[base]{NULL}}, all units are read.}
 \item{readMap}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying a read map used to read the
   cell signals.}
 \item{...}{Arguments passed to the low-level function for read units,
   e.g. \code{\link[affxparser]{readCelUnits}} or \code{\link[aroma.apd]{readApdUnits}}.}
 \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
   May also be a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns a named \code{\link[base]{list}} structure.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
