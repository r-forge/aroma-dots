%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.getSiblings.R
% 
% on Tue Nov 20 12:39:08 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getSiblings.AffymetrixCelSet}
\alias{getSiblings.AffymetrixCelSet}
\alias{AffymetrixCelSet.getSiblings}
\alias{getSiblings.AffymetrixCelSet}
\alias{getSiblings,AffymetrixCelSet-method}


\title{Gets the all data sets that refers to the same samples a this one}

\description{
 Gets the all data sets that refers to the same samples a this one.
}

\usage{\method{getSiblings}{AffymetrixCelSet}(this, notSelf=FALSE, ...)}

\arguments{
  \item{notSelf}{If \code{\link[base:logical]{TRUE}}, the current data set is not returned, just
     its siblings.}
  \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:s.AffymetrixCelSet]{*s}()} "fromFiles".}
}

\value{
 Returns a named \code{\link[base]{list}} of \code{\link{AffymetrixCelSet}} objects.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
