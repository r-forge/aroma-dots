%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  HudsonAlphaGenotypeTcgaDataFile.exportGenotypeCalls.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{exportGenotypeCalls.HudsonAlphaGenotypeTcgaDataFile}
\alias{exportGenotypeCalls.HudsonAlphaGenotypeTcgaDataFile}
\alias{HudsonAlphaGenotypeTcgaDataFile.exportGenotypeCalls}
\alias{exportGenotypeCalls.HudsonAlphaGenotypeTcgaDataFile}
\alias{exportGenotypeCalls,HudsonAlphaGenotypeTcgaDataFile-method}


\title{Exports genotype calls and confidence scores}


\description{
   Export genotype calls and corresponding confidence scores for one sample
   from a BeadStudio data file to an
    \code{\link[aroma.core]{AromaUnitGenotypeCallFile}} (genotype calls).
 }

\arguments{
  \item{dataSet}{\code{\link[base]{character}}, data set name.}
  \item{unf}{A \code{\link[aroma.core]{UnitNamesFile}}.}
  \item{...}{arguments passed to
     \code{AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(...)}.}
  \item{rootPath}{\code{\link[base]{character}}, root path for output files.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, existing output files are rewritten.}
  \item{verbose}{\code{\link{Verbose}} object.}
}


\author{
  Pierre Neuvial.
}


\details{
  A \code{\link[aroma.core]{UnitNamesFile}} is inferred from "chipType".
  Genotyping units (SNPs) are ordered according to this
  \code{\link[aroma.core]{UnitNamesFile}}.
}
\keyword{internal}
\keyword{methods}
