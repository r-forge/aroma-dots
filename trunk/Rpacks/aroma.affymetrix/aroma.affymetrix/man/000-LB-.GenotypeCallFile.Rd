%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[.GenotypeCallFile}
\alias{[.GenotypeCallFile}
\alias{GenotypeCallFile.[}
\alias{[.GenotypeCallFile}
\alias{[,GenotypeCallFile-method}


\title{Gets the genotype calls}

\description{
 Gets the genotype calls for a subset of units.
}

\usage{\method{[}{GenotypeCallFile}(this, i=NULL, drop=TRUE, ...)}

\arguments{
  \item{i}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of units to query.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{drop}{If \code{\link[base:logical]{TRUE}}, single dimensions are dropped.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{factor}} \code{\link[base]{vector}} with levels \code{-}, \code{AA}, \code{AB},
 \code{BB}, and \code{NC}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[aroma.affymetrix:readUnits.GenotypeCallFile]{*readUnits}()} is used.
  For more information see \code{\link{GenotypeCallFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
