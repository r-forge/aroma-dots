%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readUnits.GenotypeCallFile}
\alias{readUnits.GenotypeCallFile}
\alias{GenotypeCallFile.readUnits}
\alias{readUnits.GenotypeCallFile}
\alias{readUnits,GenotypeCallFile-method}


\title{Gets the genotype calls}

\description{
 Gets the genotype calls for a subset of the units.
}

\usage{\method{readUnits}{GenotypeCallFile}(this, units=NULL, ...)}

\arguments{
  \item{units}{An \code{\link[base]{vector}} specifying the subset of units to update.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{factor}} \code{\link[base]{vector}} with levels \code{-}, \code{AA}, \code{AB},
 \code{BB}, and \code{NC}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:updateUnits.GenotypeCallFile]{*updateUnits}()}.
  For more information see \code{\link{GenotypeCallFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
