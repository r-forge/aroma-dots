%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getIntensities.AffymetrixCelSet}
\alias{getIntensities.AffymetrixCelSet}
\alias{AffymetrixCelSet.getIntensities}
\alias{getIntensities.AffymetrixCelSet}
\alias{getIntensities,AffymetrixCelSet-method}


\title{Gets cell intensities from a set of cells and a set of arrays}

\description{
  Gets cell intensities from a set of cells and a set of arrays.
}

\usage{\method{getIntensities}{AffymetrixCelSet}(this, indices=NULL, verbose=FALSE, ...)}

\arguments{
 \item{arrays}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying arrays to be read.
   If \code{\link[base]{NULL}}, all arrays are read.}
 \item{indices}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying cells to be read.
   If \code{\link[base]{NULL}}, all cells are read.}
 \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
   May also be a \code{\link[R.utils]{Verbose}} object.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{numeric}} \eqn{NxK} matrix, where \eqn{N} is the number of
  cells read, and \eqn{K} is the number of arrays read.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
