%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.getAllelePairs.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSnpNames.AffymetrixCdfFile}
\alias{getSnpNames.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getSnpNames}
\alias{getSnpNames.AffymetrixCdfFile}
\alias{getSnpNames,AffymetrixCdfFile-method}


\title{Gets the names of the SNPs in the data set}

\description{
  Gets the names of the SNPs in the data set.
}

\usage{\method{getSnpNames}{AffymetrixCdfFile}(this, pattern="^SNP_A-", ...)}

\arguments{
 \item{pattern}{A regular expression to identify unit that are SNPs.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, \code{\link[aroma.affymetrix:getUnitNames.AffymetrixCdfFile]{*getUnitNames}()}.
  is used.

  For more information see \code{\link{AffymetrixCdfFile}}.
}
\keyword{internal}
\keyword{methods}
