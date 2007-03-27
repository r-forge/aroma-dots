%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectSet.extractMatrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractMatrix.ChipEffectSet}
\alias{extractMatrix.ChipEffectSet}
\alias{ChipEffectSet.extractMatrix}
\alias{extractMatrix.ChipEffectSet}
\alias{extractMatrix,ChipEffectSet-method}


\title{Extract data as a matrix for a set of arrays}

\description{
 Extract data as a matrix for a set of arrays.
}

\usage{\method{extractMatrix}{ChipEffectSet}(this, units=NULL, ..., field=c("theta", "sdTheta"), returnUgcMap=FALSE, verbose=FALSE)}

\arguments{
  \item{units}{(The subset of units to be matched.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
  \item{field}{The field to be extracted.}
  \item{returnUgcMap}{If \code{\link[base:logical]{TRUE}}, the (unit, group, cell) map is returned
    as an attribute.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an Jx(3+K) \code{\link[base]{double}} \code{\link[base]{matrix}} where J is the number of units,
 and K is the number of arrays.
 The names of the columns (after the 3rd) are the names of the arrays.
 The names of the rows are the unit indices (as indexed by the CDF).
 The rows are ordered according to \code{units} arguments.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
