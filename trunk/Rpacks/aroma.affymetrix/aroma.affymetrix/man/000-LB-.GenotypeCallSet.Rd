%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[.GenotypeCallSet}
\alias{[.GenotypeCallSet}
\alias{GenotypeCallSet.[}
\alias{[.GenotypeCallSet}
\alias{[,GenotypeCallSet-method}

\alias{GenotypeCallSet.[[}
\alias{[[.GenotypeCallSet}
\alias{[[,GenotypeCallSet-method}


\title{Gets the genotype calls}

\description{
 Gets the genotype calls for a subset of units and a subset of arrays.
}

\usage{\method{[}{GenotypeCallSet}(this, i=NULL, j=NULL, ..., drop=FALSE)}

\arguments{
  \item{i}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of units to query.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{j}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of arrays to query.
    If \code{\link[base]{NULL}}, all arrays are considered.}
  \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:readUnits.GenotypeCallSet]{*readUnits}()}.}
  \item{drop}{If \code{\link[base:logical]{TRUE}}, single dimensions are dropped.}
}

\value{
 Returns a \code{\link[base]{factor}} \code{\link[base]{matrix}} (\code{\link[base]{vector}})
 with levels \code{-}, \code{AA}, \code{AB}, \code{BB}, and \code{NC}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[aroma.affymetrix:readUnits.GenotypeCallSet]{*readUnits}()} is used.
  For more information see \code{\link{GenotypeCallSet}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
