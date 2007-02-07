%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipGenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromChipType.DChipGenomeInformation}
\alias{fromChipType.DChipGenomeInformation}
\alias{DChipGenomeInformation.fromChipType}
\alias{fromChipType.DChipGenomeInformation}
\alias{fromChipType,DChipGenomeInformation-method}


\title{Defines a DChipGenomeInformation object by chip type}

\description{
 Defines a DChipGenomeInformation object by chip type.
}

\usage{\method{fromChipType}{DChipGenomeInformation}(static, chipType, version=NULL, ...)}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{version}{An optional \code{\link[base]{character}} string specifying the version
   string, if more than one version is available.}
 \item{...}{Not used.}
}

\value{
 Returns an \code{\link{DChipGenomeInformation}} object.
 If no file was not found, an error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{DChipGenomeInformation}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
