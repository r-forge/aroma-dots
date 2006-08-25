%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFiles.AffymetrixFileSet}
\alias{fromFiles.AffymetrixFileSet}
\alias{AffymetrixFileSet.fromFiles}
\alias{fromFiles.AffymetrixFileSet}
\alias{fromFiles,AffymetrixFileSet-method}


\title{Defines an AffymetrixFileSet object by searching for Affymetrix files}

\description{
  Defines an AffymetrixFileSet object by searching for Affymetrix files.
}

\synopsis{fromFiles.AffymetrixFileSet(static, path=NULL, pattern=NULL, recursive=TRUE, fileClass="AffymetrixFile", ...)}
\usage{AffymetrixFileSet$fromFiles(path=NULL, pattern=NULL, recursive=TRUE, fileClass="AffymetrixFile", ...)}

\arguments{
 \item{path}{The path where to search for Affymetrix files.}
 \item{pattern}{The filename pattern for match files.
    If \code{\link[base]{NULL}}, filename extensions corresponding to known subclasses
    of the abstract \code{\link{AffymetrixFile}} class are search for.}
 \item{recursive}{If \code{\link[base:logical]{TRUE}}, subdirectories are search recursively,
    otherwise not.}
 \item{...}{Other arguments passed to \code{\link[base]{list.files}}().}
}

\value{
  Returns an \code{\link{AffymetrixFileSet}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
