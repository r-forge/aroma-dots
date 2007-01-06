%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnChipEffectFile.fitGlad.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitGlad.CnChipEffectFile}
\alias{fitGlad.CnChipEffectFile}
\alias{CnChipEffectFile.fitGlad}
\alias{fitGlad.CnChipEffectFile}
\alias{fitGlad,CnChipEffectFile-method}


\title{Fits the GLAD model to copy-number estimates}

\description{
 Fits the GLAD model to copy-number estimates of a certain chromosome.
}

\usage{\method{fitGlad}{CnChipEffectFile}(this, reference, chromosomes=c(1:22, "X"), units=NULL, useStddvs=TRUE, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{reference}{The \code{\link{CnChipEffectFile}} object used as the
    reference chip effects when calculating the raw (relative) copy-number
    estimates..}
  \item{chromosomes}{The chromosomes for which the model should be fitted.}
  \item{units}{(Optional) The subset of units to be matched.}
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, any in-memory cached results are ignored.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[GLAD]{profileCGH}} object returned by \code{\link[GLAD]{glad}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[GLAD]{glad}} is used.
  For more information see \code{\link{CnChipEffectFile}}.
}
\keyword{internal}
\keyword{methods}
