%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Weights.R
% 
% on Tue Jan 15 18:37:36 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{unite.MicroarrayWeights}
\alias{unite.MicroarrayWeights}
\alias{MicroarrayWeights.unite}
\alias{unite.MicroarrayWeights}
\alias{unite,MicroarrayWeights-method}


\title{Coalesce weights of lower order into higher order weights}

\synopsis{unite.MicroarrayWeights(static, weights, ...)}
\usage{MicroarrayWeights$unite(weights, ...)}

\description{
 Coalesce weights of lower order into higher order weights.

 Each new weight is the square root of the average sum of squares of
 lower-order weights.
}

\arguments{
 \item{weights}{A NxK \code{\link[base]{matrix}} of weights to be united row by row.}
}

\value{
  Returns a weight \code{\link[base]{vector}} of length N.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MicroarrayWeights}}..
}
\keyword{internal}
\keyword{methods}
