%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectSet.fitCrlmm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitCrlmm.SnpChipEffectSet}
\alias{fitCrlmm.SnpChipEffectSet}
\alias{SnpChipEffectSet.fitCrlmm}
\alias{fitCrlmm.SnpChipEffectSet}
\alias{fitCrlmm,SnpChipEffectSet-method}


\title{Fits chip effects to the CRLMM model}

\description{
 Fits chip effects to the CRLMM model.
}

\usage{\method{fitCrlmm}{SnpChipEffectSet}(this, minLLRforCalls=c(AA = 50, AB = 40, BB = 50), recalibrate=TRUE, transform=c("log", "asinh"), ..., verbose=FALSE)}

\arguments{
  \item{minLLRforCalls}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length three.}
  \item{recalibrate}{If \code{\link[base:logical]{TRUE}}, a second round of the fitting is done.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[oligo]{SnpCallSetPlus-class}} object.
}

\details{
  Simple benchmarking on shadowfax:
  For 90 CEPH Xba chips it takes ~90 minutes, i.e. 60 s/array.
}

\author{Henrik Bengtsson \url{http://www.braju.com/R/}}

\seealso{
  \code{\link[oligo]{crlmm}}.
  For more information see \code{\link{SnpChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
