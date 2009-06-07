%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.ARGS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{validateArgumentWeights.MicroarrayData}
\alias{validateArgumentWeights.MicroarrayData}
\alias{MicroarrayData.validateArgumentWeights}
\alias{validateArgumentWeights.MicroarrayData}
\alias{validateArgumentWeights,MicroarrayData-method}


\title{Validates the argument 'weights'}

\description{
 Validates the argument 'weights'.
}

\usage{\method{validateArgumentWeights}{MicroarrayData}(this, weights, slides=NULL, zeroOneOnly=FALSE, ...)}

\arguments{
  \item{weights}{Either a \code{\link[base]{numeric}}, \code{\link[base]{logical}} \code{\link[base]{vector}}, or \code{\link[base]{NULL}}.
    If a \code{\link[base]{logical}}, it is coersed to a zero-one \code{\link[base]{double}} \code{\link[base]{vector}}.
    If a \code{\link[base]{numeric}}, each element must be in [0,1], otherwise an \code{\link[R.oo]{Exception}}
    is thrown.
    If \code{\link[base]{NULL}}, then \code{\link[base]{NULL}} is returned.
  }
  \item{zeroOneOnly}{If \code{\link[base:logical]{TRUE}}, only zero-one weights are allows.
    All non-zero weights are set to one.
  }
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}} of weights of equal length as the input \code{\link[base]{vector}}.
 If \code{\link[base]{NULL}} was given, then \code{\link[base]{NULL}} is returned.

 If the argument was invalid an \code{\link[R.oo]{Exception}} is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
