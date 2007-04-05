%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnChipEffectSet.importFromDChip.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{importFromDChip.CnChipEffectSet}
\alias{importFromDChip.CnChipEffectSet}
\alias{CnChipEffectSet.importFromDChip}
\alias{importFromDChip.CnChipEffectSet}
\alias{importFromDChip,CnChipEffectSet-method}


\title{Imports copy-number estimates from a dChip result file}

\description{
 Imports copy-number estimates from a dChip result file.
 Currently only total copy-number estimates can be imported, that is
 if dChip fitted the PLM with allele A and allele B combined.
}

\synopsis{importFromDChip.CnChipEffectSet(static, filename, path=NULL, combineAlleles=TRUE, cdf, ..., skip=TRUE, verbose=FALSE)}
\usage{CnChipEffectSet$importFromDChip(filename, path=NULL, combineAlleles=TRUE, cdf, ..., skip=TRUE, verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the dChip result file.}
  \item{path}{An optional path to the file.}
  \item{combineAlleles}{If \code{\link[base:logical]{TRUE}}, .}
  \item{cdf}{An \code{\link{AffymetrixCdfFile}} object.}
  \item{...}{Not used.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, already imported arrays will be skipped.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{CnChipEffectSet}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CnChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
