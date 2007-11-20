%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% on Tue Nov 20 12:38:56 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findByChipType.AffymetrixCdfFile}
\alias{findByChipType.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.findByChipType}
\alias{findByChipType.AffymetrixCdfFile}
\alias{findByChipType,AffymetrixCdfFile-method}


\title{Locates a CDF file from its chip type}

\description{
 Locates a CDF file from its chip type.
}

\synopsis{findByChipType.AffymetrixCdfFile(static, chipType, pattern=NULL, ..., .useAffxparser=TRUE)}
\usage{AffymetrixCdfFile$findByChipType(chipType, pattern=NULL, ..., .useAffxparser=TRUE)}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{...}{Not used.}
 \item{.useAffxparser}{If \code{\link[base:logical]{TRUE}}, \code{\link[affxparser]{findCdf}} is used if
   the CDF could not be located.}
}

\value{
 Returns a pathname as a \code{\link[base]{character}} string to the first CDF file found.
 If non CDF with requested chip type was found, \code{\link[base]{NULL}} is returned.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\details{
}

\seealso{
  \code{\link[aroma.affymetrix:byChipType.AffymetrixCdfFile]{*byChipType}()}.
  For more information see \code{\link{AffymetrixCdfFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
