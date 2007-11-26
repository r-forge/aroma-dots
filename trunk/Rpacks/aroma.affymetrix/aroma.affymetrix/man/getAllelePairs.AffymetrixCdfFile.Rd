%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.getAllelePairs.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getAllelePairs.AffymetrixCdfFile}
\alias{getAllelePairs.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getAllelePairs}
\alias{getAllelePairs.AffymetrixCdfFile}
\alias{getAllelePairs,AffymetrixCdfFile-method}


\title{Gets the allele pairs for each unit}

\description{
  Gets the allele pairs for each unit.
}

\usage{\method{getAllelePairs}{AffymetrixCdfFile}(this, units=NULL, ...)}

\arguments{
 \item{units}{Subset of units to be queried.  If \code{\link[base]{NULL}}, all units are used.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{vector}} of \code{\link[base]{factor}}s.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getAllelePairUnitSets.AffymetrixCdfFile]{*getAllelePairUnitSets}()}.
  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
