%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSetTuple.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{extract.AffymetrixCelSetTuple}
\alias{extract.AffymetrixCelSetTuple}
\alias{AffymetrixCelSetTuple.extract}
\alias{extract.AffymetrixCelSetTuple}
\alias{extract,AffymetrixCelSetTuple-method}


\title{Extracts a subset AffymetrixCelSetTuple}

\description{
 Extracts a subset AffymetrixCelSetTuple.
}

\usage{\method{extract}{AffymetrixCelSetTuple}(this, arrays, ..., verbose=FALSE)}

\arguments{
  \item{arrays}{A \code{\link[base]{vector}} of arrays to be extracted.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link{AffymetrixCelSetTuple}}.
}

\details{
  If no matching arrays are available for a given chip type, that chip type
  is excluded in the returned set tuple.  This guarantees that a set tuple
  always contains existing arrays.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSetTuple}}.
}
\keyword{internal}
\keyword{methods}
