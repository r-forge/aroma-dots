%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[< -.GenotypeCallFile}
\alias{[<-.GenotypeCallFile}
\alias{GenotypeCallFile.[<-}
\alias{[<-.GenotypeCallFile}
\alias{[<-,GenotypeCallFile-method}


\title{Sets the genotype calls}

\description{
 Sets the genotype calls for a subset of the units.
}

\usage{\method{[}{GenotypeCallFile}(this, i=NULL) <- value}

\arguments{
  \item{i}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the subset of units to update.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{value}{A \code{\link[base]{vector}} of the same length as \code{i} with genotype
    calls.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[aroma.affymetrix:updateUnits.GenotypeCallFile]{*updateUnits}()} is used.
  For more information see \code{\link{GenotypeCallFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
