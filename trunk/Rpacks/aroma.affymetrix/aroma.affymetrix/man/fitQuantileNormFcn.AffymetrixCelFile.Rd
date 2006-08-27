%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.normalizeQuantiles.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fitQuantileNormFcn.AffymetrixCelFile}
\alias{fitQuantileNormFcn.AffymetrixCelFile}
\alias{AffymetrixCelFile.fitQuantileNormFcn}
\alias{fitQuantileNormFcn.AffymetrixCelFile}
\alias{fitQuantileNormFcn,AffymetrixCelFile-method}


\title{Fits quantile normalization functions for the arrays in the data set}

\description{
 Fits quantile normalization functions for the arrays in the data set.
}

\usage{\method{fitQuantileNormFcn}{AffymetrixCelFile}(this, yTarget, subset=NULL, ..., controlParams=list(spar = NULL, nknots = 1024), verbose=FALSE)}

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
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
