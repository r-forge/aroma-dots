%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromChipType.AffymetrixCdfFile}
\alias{fromChipType.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.fromChipType}
\alias{fromChipType.AffymetrixCdfFile}
\alias{fromChipType,AffymetrixCdfFile-method}


\title{Defines an AffymetrixCdfFile object by chip type}

\description{
 Defines an AffymetrixCdfFile object by chip type utilizing the \code{\link[affxparser]{findCdf}} method.
}

\synopsis{fromChipType.AffymetrixCdfFile(static, chipType, ...)}
\usage{AffymetrixCdfFile$fromChipType(chipType, ...)}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{...}{Not used.}
}

\value{
 Returns an \code{\link{AffymetrixCdfFile}} object.
 If a matching CDF file was not found, an error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:fromFile.AffymetrixCdfFile]{*fromFile}()}.
  For more information see \code{\link{AffymetrixCdfFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
