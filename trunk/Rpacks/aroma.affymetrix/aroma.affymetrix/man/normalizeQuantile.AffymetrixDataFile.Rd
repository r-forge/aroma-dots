%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantile.AffymetrixDataFile}
\alias{normalizeQuantile.AffymetrixDataFile}
\alias{AffymetrixDataFile.normalizeQuantile}
\alias{normalizeQuantile.AffymetrixDataFile}
\alias{normalizeQuantile,AffymetrixDataFile-method}


\title{Normalizes the sample to a target empirical distribution}

\description{
 Normalizes the sample to a target empirical distribution.
}

\usage{\method{normalizeQuantile}{AffymetrixDataFile}(this, outPath=file.path("norm", getChipType(this)), xTarget, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, format=getOption("aroma.affymetrix/format"), verbose=FALSE)}

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
  \item{format}{The output file format.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the normalized \code{\link{AffymetrixDataFile}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.light]{normalizeQuantile.numeric}}
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
