%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GladModel.fitOne.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitOne.GladModel}
\alias{fitOne.GladModel}
\alias{GladModel.fitOne}
\alias{fitOne.GladModel}
\alias{fitOne,GladModel-method}


\title{Fits the GLAD model for one chromosome in one sample}

\description{
 Fits the GLAD model for one chromosome in one sample.
}

\usage{\method{fitOne}{GladModel}(this, ceList, refList, chromosome, units=NULL, useStddvs=TRUE, ..., maxNAFraction=1/8, force=FALSE, verbose=FALSE)}

\arguments{
  \item{ceList}{A \code{\link[base]{list}} of \code{\link{ChipEffectFile}} objects.}
  \item{refList}{A \code{\link[base]{list}} of \code{\link{ChipEffectFile}} objects.}
  \item{chromosome}{The chromosome for which the model should be fitted.}
  \item{units}{(Optional) The subset of units to be matched.}
  \item{useStdvs}{If \code{\link[base:logical]{TRUE}}, standard deviations estimates for the
     chip effects are passed along and returned.  However, note that
     GLAD is not making use of these estimates.}
  \item{...}{Not used.}
  \item{maxNAFraction}{A \code{\link[base]{numeric}} in [0,1] specifying the maximum fraction
     of non-finite values allowed.
     If more are detected, this is interpreted as something has gone wrong
     in the preprocessing and an error is thrown.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, any in-memory cached results are ignored.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[GLAD]{profileCGH}} object returned by \code{\link[GLAD]{glad}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[GLAD]{glad}} is used.
  For more information see \code{\link{GladModel}}.
}
\keyword{internal}
\keyword{methods}
