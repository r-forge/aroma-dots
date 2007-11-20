%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% on Tue Nov 20 12:38:56 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitSizes.AffymetrixCdfFile}
\alias{getUnitSizes.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getUnitSizes}
\alias{getUnitSizes.AffymetrixCdfFile}
\alias{getUnitSizes,AffymetrixCdfFile-method}


\title{Gets the number of groups in each unit}

\description{
 Gets the number of groups in each unit.
}

\usage{\method{getUnitSizes}{AffymetrixCdfFile}(this, units=NULL, ...)}

\arguments{
  \item{units}{The units of interest. If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{vector}} of \code{\link[base]{integer}}s.
}

\details{
  Once read from file, this information is cached in memory for efficiency.
  The cache can be cleared by calling \code{gc(cdf)}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
