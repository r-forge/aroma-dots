%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  exportGenotypeCallsAndConfidenceScores.BroadBirdseedGenotypeTcgaDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{exportGenotypeCallsAndConfidenceScores.BroadBirdseedGenotypeTcgaDataFile}
\alias{exportGenotypeCallsAndConfidenceScores.BroadBirdseedGenotypeTcgaDataFile}
\alias{BroadBirdseedGenotypeTcgaDataFile.exportGenotypeCallsAndConfidenceScores}
\alias{exportGenotypeCallsAndConfidenceScores.BroadBirdseedGenotypeTcgaDataFile}
\alias{exportGenotypeCallsAndConfidenceScores,BroadBirdseedGenotypeTcgaDataFile-method}


\title{Exports genotype calls and confidence scores}


\description{
   Export genotype calls and corresponding confidence scores for one sample
   from a Birdseed data file to an
    \code{\link[aroma.core]{AromaUnitGenotypeCallFile}} (genotype calls) and
    an \code{\link[aroma.core]{AromaUnitSignalBinaryFile}} (confidence scores).
 }

\arguments{
  \item{dataSet}{\code{\link[base]{character}}, data set name.}
  \item{chipType}{\code{\link[base]{character}}, chip type.}
  \item{...}{arguments passed to
     \code{AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(...)}.}
  \item{rootPath}{\code{\link[base]{character}}, root path for output files.}
  \item{verbose}{\code{\link{Verbose}} object.}
}


\author{
  Pierre Neuvial.
}


\details{
  A \code{\link[aroma.core]{UnitNamesFile}} is inferred from "chipType" using
  \code{AffymetrixCdfFile$byChipType()}.
  Genotyping units (SNPs) are ordered according to this
  \code{\link[aroma.core]{UnitNamesFile}}.
  Confidence scores are stored as one minus the input score as the
  latter is a p-value.
}

\note{
  The implementation requires \pkg{aroma.affymetrix} only for the
  inference of a \code{\link[aroma.core]{UnitNamesFile}}.  However it is tied
  to Affymetrix data anyway because Birdseed only supports Affymetrix data.
}
\keyword{internal}
\keyword{methods}
