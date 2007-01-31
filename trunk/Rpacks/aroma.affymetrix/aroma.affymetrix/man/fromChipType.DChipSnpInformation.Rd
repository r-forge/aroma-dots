%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipSnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromChipType.DChipSnpInformation}
\alias{fromChipType.DChipSnpInformation}
\alias{DChipSnpInformation.fromChipType}
\alias{fromChipType.DChipSnpInformation}
\alias{fromChipType,DChipSnpInformation-method}


\title{Defines a DChipSnpInformation object by chip type}

\description{
 Defines a DChipSnpInformation object by chip type.
}

\usage{\method{fromChipType}{DChipSnpInformation}(static, chipType, rootPath="annotations", version=NULL, pattern=NULL, ...)}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{rootPath}{A \code{\link[base]{character}} string specifying the root path, i.e.
   the annotation directory.}
 \item{version}{An optional \code{\link[base]{character}} string specifying the version
   string, if more than one version is available.}
 \item{pattern}{An optional filename pattern used to locate the
   dChip genome file.  If \code{\link[base]{NULL}}, a default pattern is used.}
 \item{...}{Not used.}
}

\value{
 Returns an \code{\link{DChipSnpInformation}} object.
 If no file was not found, an error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{DChipSnpInformation}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
