%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{updateUnits.GenotypeCallFile}
\alias{updateUnits.GenotypeCallFile}
\alias{GenotypeCallFile.updateUnits}
\alias{updateUnits.GenotypeCallFile}
\alias{updateUnits,GenotypeCallFile-method}


\title{Updates the genotype calls}

\description{
 Updates the genotype calls for a subset of the units.
}

\usage{\method{updateUnits}{GenotypeCallFile}(this, units=NULL, calls, ...)}

\arguments{
  \item{units}{An \code{\link[base]{vector}} specifying the subset of units to update.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{calls}{A \code{\link[base]{vector}} of the same length as \code{units} with
    genotype calls.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:readUnits.GenotypeCallFile]{*readUnits}()}.
  For more information see \code{\link{GenotypeCallFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
