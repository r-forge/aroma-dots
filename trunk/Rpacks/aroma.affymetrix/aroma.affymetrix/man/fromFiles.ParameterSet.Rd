%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFiles.ParameterSet}
\alias{fromFiles.ParameterSet}
\alias{ParameterSet.fromFiles}
\alias{fromFiles.ParameterSet}
\alias{fromFiles,ParameterSet-method}


\title{Defines a ParameterSet object by searching for files}

\description{
  Defines a ParameterSet object by searching for files.
}

\synopsis{fromFiles.ParameterSet(static, path="estimates/", pattern=NULL, recursive=TRUE, ...)}
\usage{ParameterSet$fromFiles(path="estimates/", pattern=NULL, recursive=TRUE, ...)}

\arguments{
 \item{path}{The path where to search for files.}
 \item{pattern}{The filename pattern for match files.
    If \code{\link[base]{NULL}}, filename extensions corresponding to known subclasses
    of the abstract \code{\link{ParameterFile}} class are search for.}
 \item{recursive}{If \code{\link[base:logical]{TRUE}}, subdirectories are search recursively,
    otherwise not.}
 \item{...}{Other arguments passed to \code{\link[base]{list.files}}().}
}

\value{
  Returns an \code{\link{ParameterSet}} object.
}

\author{Henrik Bengtsson \url{http://www.braju.com/R/}}

\seealso{
  For more information see \code{\link{ParameterSet}}.
}
\keyword{internal}
\keyword{methods}
