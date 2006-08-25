%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fitQuantileNormFcn.AffymetrixDataFile}
\alias{fitQuantileNormFcn.AffymetrixDataFile}
\alias{AffymetrixDataFile.fitQuantileNormFcn}
\alias{fitQuantileNormFcn.AffymetrixDataFile}
\alias{fitQuantileNormFcn,AffymetrixDataFile-method}


\title{Fits quantile normalization functions for the arrays in the data set}

\description{
 Fits quantile normalization functions for the arrays in the data set.
}

\usage{\method{fitQuantileNormFcn}{AffymetrixDataFile}(this, yTarget, subset=NULL, ..., controlParams=list(spar = NULL, nknots = 1024), verbose=FALSE)}

\arguments{
  \item{yTarget}{The target probe signals.}
  \item{subset}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying a subset of probe
     indices used to fit the normalization function.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a normalization \code{\link[base]{function}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataFile}}.
}
\keyword{internal}
\keyword{methods}
