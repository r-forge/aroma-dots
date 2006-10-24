%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFileType.AffymetrixFile}
\alias{getFileType.AffymetrixFile}
\alias{AffymetrixFile.getFileType}
\alias{getFileType.AffymetrixFile}
\alias{getFileType,AffymetrixFile-method}


\title{Gets the file type of a file}

\description{
  Gets the file type of a file.
}

\usage{\method{getFileType}{AffymetrixFile}(this, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}} in lower case letters.
}

\details{
  By default, this methods returns the filename extension, but subclasses
  may override this.
}

\author{Henrik Bengtsson \url{http://www.braju.com/R/}}

\seealso{
  For more information see \code{\link{AffymetrixFile}}.
}
\keyword{internal}
\keyword{methods}
