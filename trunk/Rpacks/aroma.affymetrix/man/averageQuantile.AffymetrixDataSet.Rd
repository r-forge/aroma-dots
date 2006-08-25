%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{averageQuantile.AffymetrixDataSet}
\alias{averageQuantile.AffymetrixDataSet}
\alias{AffymetrixDataSet.averageQuantile}
\alias{averageQuantile.AffymetrixDataSet}
\alias{averageQuantile,AffymetrixDataSet-method}


\title{Gets the average empirical distribution across all samples}

\description{
 Gets the average empirical distribution across all samples.
}

\usage{\method{averageQuantile}{AffymetrixDataSet}(this, probes=NULL, ..., verbose=FALSE)}

\arguments{
  \item{probes}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying what subset of
     probes to be used to calculate the empirical distribution.
     If \code{\link[base]{NULL}}, all probes are used.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.light]{averageQuantile.list}}
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
