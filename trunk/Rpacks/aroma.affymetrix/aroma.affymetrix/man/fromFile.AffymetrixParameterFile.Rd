%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixParameterFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFile.AffymetrixParameterFile}
\alias{fromFile.AffymetrixParameterFile}
\alias{AffymetrixParameterFile.fromFile}
\alias{fromFile.AffymetrixParameterFile}
\alias{fromFile,AffymetrixParameterFile-method}


\title{Defines an AffymetrixParameterFile object from a file}

\description{
 Defines an AffymetrixParameterFile object from a file.
}

\synopsis{fromFile.AffymetrixParameterFile(static, filename, path=NULL, ..., .checkArgs=TRUE)}
\usage{AffymetrixParameterFile$fromFile(filename, path=NULL, ..., .checkArgs=TRUE)}

\arguments{
 \item{filename}{The filename of to the file.}
 \item{path}{The path to the file.}
 \item{...}{Not used.}
}

\value{
 Returns an \code{\link{AffymetrixParameterFile}} object.
 If the file is not found or if it is of the wrong file format, an
 error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixParameterFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
