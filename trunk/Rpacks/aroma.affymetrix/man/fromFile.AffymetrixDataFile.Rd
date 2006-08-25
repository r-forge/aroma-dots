%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.IO.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fromFile.AffymetrixDataFile}
\alias{fromFile.AffymetrixDataFile}
\alias{AffymetrixDataFile.fromFile}
\alias{fromFile.AffymetrixDataFile}
\alias{fromFile,AffymetrixDataFile-method}


\title{Defines an AffymetrixDataFile object from a file}

\description{
 Defines an AffymetrixDataFile object from a file by calling the same static method of all known subclasses
 "bottom up" until the data file is "accepted".
}

\synopsis{fromFile.AffymetrixDataFile(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE)}
\usage{AffymetrixDataFile$fromFile(filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE)}

\arguments{
  \item{filename}{The filename of to the file.}
  \item{path}{The path to the file.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.checkArgs}{Internal.}
}

\value{
 Returns an instance of a subclass of \code{\link{AffymetrixDataFile}}.
 If the data file is not recognized by any of the subclasses, an
 error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
