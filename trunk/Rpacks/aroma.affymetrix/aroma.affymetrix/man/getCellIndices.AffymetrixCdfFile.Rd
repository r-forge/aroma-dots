%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% on Sat Jan 06 01:46:38 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCellIndices.AffymetrixCdfFile}
\alias{getCellIndices.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getCellIndices}
\alias{getCellIndices.AffymetrixCdfFile}
\alias{getCellIndices,AffymetrixCdfFile-method}


\title{Gets the cell indices unit by unit}

\description{
 Gets the cell indices unit by unit of all or a subset of the units.
}

\usage{\method{getCellIndices}{AffymetrixCdfFile}(this, units=NULL, ..., verbose=FALSE)}

\arguments{
  \item{units}{The units of interest. If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Additional arguments passed to
     \code{\link[affxparser]{readCdfCellIndices}}.}
}

\value{
 Returns the \code{\link[base]{list}} structure returned by
 \code{\link[affxparser]{readCdfCellIndices}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  See \code{\link[aroma.affymetrix:setRestructor.AffymetrixCdfFile]{*setRestructor}()} to set a default re-constructor for
  the returned CDF structure.
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
