%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectSet.calculateBaseline.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{calculateBaseline.ChipEffectSet}
\alias{calculateBaseline.ChipEffectSet}
\alias{ChipEffectSet.calculateBaseline}
\alias{calculateBaseline.ChipEffectSet}
\alias{calculateBaseline,ChipEffectSet-method}


\title{Estimates the baseline signal chromosome by chromosome}

\description{
 Estimates the baseline signal chromosome by chromosome.
}

\usage{\method{calculateBaseline}{ChipEffectSet}(this, chromosomes=NULL, ploidy=2, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{chromosomes}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying for which chromsosomes
    the baseline should be estimated.
    If \code{\link[base]{NULL}}, all chromosomes are considered.}
  \item{ploidy}{An \code{\link[base]{integer}} specifying the ploidy that the baseline
    should have.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
