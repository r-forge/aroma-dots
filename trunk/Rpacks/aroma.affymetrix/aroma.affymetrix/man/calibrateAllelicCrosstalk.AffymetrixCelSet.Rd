%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.calibrateAllelicCrosstalk.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{calibrateAllelicCrosstalk.AffymetrixCelSet}
\alias{calibrateAllelicCrosstalk.AffymetrixCelSet}
\alias{AffymetrixCelSet.calibrateAllelicCrosstalk}
\alias{calibrateAllelicCrosstalk.AffymetrixCelSet}
\alias{calibrateAllelicCrosstalk,AffymetrixCelSet-method}


\title{Calibrates probepair signals for crosstalk between allele A and allele B}

\description{
 Calibrates probepair signals for crosstalk between allele A and allele B.
}

\usage{\method{calibrateAllelicCrosstalk}{AffymetrixCelSet}(this, path=NULL, name="calibAllelicCT", ..., verbose=FALSE)}

\arguments{
  \item{outPath}{The path where to save the calibrated data files.}
  \item{...}{Additional arguments passed to
    \code{\link[aroma.affymetrix:calibrateAllelicCrosstalk.AffymetrixCelSet]{*calibrateAllelicCrosstalk}()}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the calibrated \code{\link{AffymetrixCelSet}} object.
}

\details{
 The calibration is done for each array separately.
}

\section{Benchmarking}{
  On an IBM Thinkpad A31 with 1.8GHz and 1GB RAM:
  \itemize{
   \item{Mapping50K\_Xba240}{58960 SNPs: ~55-60s/array.}
  }
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[sfit]{cfit}}.
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
