%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UgpGenomeInformation.R
% 
% on Tue Nov 20 12:41:27 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromChipType.UgpGenomeInformation}
\alias{fromChipType.UgpGenomeInformation}
\alias{UgpGenomeInformation.fromChipType}
\alias{fromChipType.UgpGenomeInformation}
\alias{fromChipType,UgpGenomeInformation-method}


\title{Defines a UgpGenomeInformation object by chip type}

\description{
 Defines a UgpGenomeInformation object by chip type.
}

\synopsis{fromChipType.UgpGenomeInformation(static, chipType, ..., verbose=FALSE)}
\usage{UgpGenomeInformation$fromChipType(chipType, ..., verbose=FALSE)}

\arguments{
 \item{chipType}{A \code{\link[base]{character}} string.}
 \item{version}{An optional \code{\link[base]{character}} string specifying the version
   string, if more than one version is available.}
 \item{...}{Not used.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{UgpGenomeInformation}} object.
 If no file was not found, an error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{UgpGenomeInformation}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
