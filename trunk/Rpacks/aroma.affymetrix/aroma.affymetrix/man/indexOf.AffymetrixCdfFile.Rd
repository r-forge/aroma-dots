%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% on Tue Nov 20 12:38:56 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{indexOf.AffymetrixCdfFile}
\alias{indexOf.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.indexOf}
\alias{indexOf.AffymetrixCdfFile}
\alias{indexOf,AffymetrixCdfFile-method}


\title{Gets the indices of units by their names}

\description{
 Gets the indices of units by their names.
}

\usage{\method{indexOf}{AffymetrixCdfFile}(this, pattern=NULL, names=NULL, ...)}

\arguments{
  \item{pattern}{A pattern to be used for identifying unit names of
     interest.  If \code{\link[base]{NULL}}, no regular expression matching is done.}
  \item{names}{Names to be match exactly to the unit names.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{vector}} of \code{\link[base]{integer}}s in [1,N] where N is the number of units
 in this CDF structure.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getUnitNames.AffymetrixCdfFile]{*getUnitNames}()}.
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
