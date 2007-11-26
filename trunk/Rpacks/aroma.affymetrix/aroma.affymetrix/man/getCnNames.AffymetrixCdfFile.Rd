%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.getAllelePairs.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCnNames.AffymetrixCdfFile}
\alias{getCnNames.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getCnNames}
\alias{getCnNames.AffymetrixCdfFile}
\alias{getCnNames,AffymetrixCdfFile-method}


\title{Gets the names of the CN units}

\description{
  Gets the names of the CN units.
}

\usage{\method{getCnNames}{AffymetrixCdfFile}(this, pattern="CN_", ...)}

\arguments{
 \item{pattern}{A regular expression to identify unit that are CN units.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getSnpNames.AffymetrixCdfFile]{*getSnpNames}()}.
  Internally, \code{\link[aroma.affymetrix:getUnitNames.AffymetrixCdfFile]{*getUnitNames}()}.
  is used.

  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
