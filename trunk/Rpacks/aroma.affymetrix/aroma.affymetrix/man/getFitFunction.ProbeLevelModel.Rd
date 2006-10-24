%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFitFunction.ProbeLevelModel}
\alias{getFitFunction.ProbeLevelModel}
\alias{ProbeLevelModel.getFitFunction}
\alias{getFitFunction.ProbeLevelModel}
\alias{getFitFunction,ProbeLevelModel-method}


\title{Static method to get the low-level function that fits the PLM}

\description{
 Static method to get the low-level function that fits the PLM.
 Any subclass model must provide this method, which should return
 a \code{\link[base]{function}} that accepts an IxJ \code{\link[base]{matrix}}.
}

\synopsis{getFitFunction.ProbeLevelModel(...)}
\usage{ProbeLevelModel$getFitFunction()}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{function}}.
}

\author{Henrik Bengtsson \url{http://www.braju.com/R/}}

\seealso{
  For more information see \code{\link{ProbeLevelModel}}.
}
\keyword{internal}
\keyword{methods}
