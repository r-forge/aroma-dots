%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSetTuple.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getChipTypes.AffymetrixCelSetTuple}
\alias{getChipTypes.AffymetrixCelSetTuple}
\alias{AffymetrixCelSetTuple.getChipTypes}
\alias{getChipTypes.AffymetrixCelSetTuple}
\alias{getChipTypes,AffymetrixCelSetTuple-method}


\title{Gets a label for all chip types merged}

\description{
 Gets a label for all chip types merged.
}

\usage{\method{getChipTypes}{AffymetrixCelSetTuple}(this, merge=FALSE, collapse="+", ...)}

\arguments{
  \item{merge}{If \code{\link[base:logical]{TRUE}}, the chip types are merged into a single string
    compressed such that only non-common parts are replicated.}
  \item{collapse}{The string used to merge chip types.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{character}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSetTuple}}.
}
\keyword{internal}
\keyword{methods}
