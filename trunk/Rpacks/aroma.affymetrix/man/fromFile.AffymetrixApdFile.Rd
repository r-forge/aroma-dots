%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixApdFile.IO.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fromFile.AffymetrixApdFile}
\alias{fromFile.AffymetrixApdFile}
\alias{AffymetrixApdFile.fromFile}
\alias{fromFile.AffymetrixApdFile}
\alias{fromFile,AffymetrixApdFile-method}


\title{Defines an AffymetrixApdFile object from a file}

\description{
 Defines an AffymetrixApdFile object from a file.
}

\synopsis{fromFile.AffymetrixApdFile(static, filename, path=NULL, ..., .checkArgs=TRUE)}
\usage{AffymetrixApdFile$fromFile(filename, path=NULL, ..., .checkArgs=TRUE)}

\arguments{
 \item{filename}{The filename of to the file.}
 \item{path}{The path to the file.}
 \item{...}{Not used.}
}

\value{
 Returns an \code{\link{AffymetrixApdFile}} object.
 If the file is not found or if it is of the wrong file format, an
 error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixApdFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
