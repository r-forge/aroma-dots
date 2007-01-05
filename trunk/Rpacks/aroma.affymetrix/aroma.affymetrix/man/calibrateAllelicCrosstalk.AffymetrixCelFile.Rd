%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.calibrateAllelicCrosstalk.R
% 
% on Sat Jan 06 01:46:41 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{calibrateAllelicCrosstalk.AffymetrixCelFile}
\alias{calibrateAllelicCrosstalk.AffymetrixCelFile}
\alias{AffymetrixCelFile.calibrateAllelicCrosstalk}
\alias{calibrateAllelicCrosstalk.AffymetrixCelFile}
\alias{calibrateAllelicCrosstalk,AffymetrixCelFile-method}


\title{Calibrates probepair signals for allele A allele B crosstalk}

\description{
 Calibrates probepair signals for allele A allele B crosstalk.
}

\usage{\method{calibrateAllelicCrosstalk}{AffymetrixCelFile}(this, path, alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, targetAvg=NULL, ..., setsOfProbes=NULL, overwrite=FALSE, skip=!overwrite, verbose=FALSE)}

\arguments{
  \item{path}{The path where to save the calibrated data files.}
  \item{alpha}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} in [0,1] of the sequence of alpha
    parameter in the cone fit.}
  \item{q}{A \code{\link[base]{numeric}} in [0,100] of the q parameter in the cone fit.}
  \item{Q}{A \code{\link[base]{numeric}} in [0,100] of the q parameter in the cone fit.}
  \item{targetAvg}{A \code{\link[base]{numeric}} specifying the average signal for all
    probes of allele A and allele B, respectively.}
  \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:getData.AffymetrixCelFile]{*getData}()}.}
  \item{setsOfProbes}{A named \code{\link[base]{list}} where each item specifies the probe indices
    for each of the four groups; forward and reverse strands on allele A
    and allele B.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already normalized arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{format}{The output file format.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the calibrated \code{\link{AffymetrixCelFile}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link{calibrateAllelicCrosstalk.AffymetrixCelSet}}().
  \code{\link[sfit]{cfit}}.
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
