%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnChipEffectSet.fitGlad.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitGlad.CnChipEffectSet}
\alias{fitGlad.CnChipEffectSet}
\alias{CnChipEffectSet.fitGlad}
\alias{fitGlad.CnChipEffectSet}
\alias{fitGlad,CnChipEffectSet-method}


\title{Fits the GLAD model to copy-number estimates}

\description{
 Fits the GLAD model to copy-number estimates of all files in the data set for a certain chromosome.
}

\usage{\method{fitGlad}{CnChipEffectSet}(this, reference=NULL, arrays=1:nbrOfArrays(this), chromosomes=getChromosomes(this), ..., verbose=FALSE)}

\arguments{
  \item{...}{Arguments passed to \code{fitGlad()} of
    \code{\link{CnChipEffectFile}}.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, any in-memory cached results are ignored.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns...
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CnChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
