%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.normalizeQuantiles.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantiles.AffymetrixCelFile}
\alias{normalizeQuantiles.AffymetrixCelFile}
\alias{AffymetrixCelFile.normalizeQuantiles}
\alias{normalizeQuantiles.AffymetrixCelFile}
\alias{normalizeQuantiles,AffymetrixCelFile-method}


\title{Normalizes the probe intensities to a target empirical distribution}

\description{
 Normalizes the probe intensities to a target empirical distribution.
}

\usage{\method{normalizeQuantiles}{AffymetrixCelFile}(this, outPath=file.path("normQuantile", getChipType(this)), xTarget, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, verbose=FALSE)}

\arguments{
  \item{outPath}{The path where to save the normalized data files.}
  \item{xTarget}{A \code{\link[base]{numeric}} \code{\link[base]{vector}}.  The empirical distribution
    to which all arrays should be normalized to.}
  \item{subsetToUpdate}{The probes to be updated.
    If \code{\link[base]{NULL}}, all probes are updated.}
  \item{typesToUpdate}{Types of probes to be updated.}
  \item{...}{Additional arguments passed to \code{normalizeQuantile()}.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already normalized arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the normalized \code{\link{AffymetrixCelFile}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.light]{normalizeQuantile.numeric}}
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
