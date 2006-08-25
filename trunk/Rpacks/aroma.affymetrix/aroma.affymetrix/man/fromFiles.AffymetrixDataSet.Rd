%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.IO.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fromFiles.AffymetrixDataSet}
\alias{fromFiles.AffymetrixDataSet}
\alias{AffymetrixDataSet.fromFiles}
\alias{fromFiles.AffymetrixDataSet}
\alias{fromFiles,AffymetrixDataSet-method}


\title{Defines an AffymetrixDataSet object by searching for Affymetrix files}

\description{
  Defines an AffymetrixDataSet object by searching for Affymetrix files.
}

\synopsis{fromFiles.AffymetrixDataSet(static, path="cel/", pattern=NULL, recursive=TRUE, ...)}
\usage{AffymetrixDataSet$fromFiles(path="cel/", pattern=NULL, recursive=TRUE, ...)}

\arguments{
 \item{path}{The path where to search for Affymetrix files.}
 \item{pattern}{The filename pattern for match files.
    If \code{\link[base]{NULL}}, filename extensions corresponding to known subclasses
    of the abstract \code{\link{AffymetrixDataFile}} class are search for.}
 \item{recursive}{If \code{\link[base:logical]{TRUE}}, subdirectories are search recursively,
    otherwise not.}
 \item{...}{Other arguments passed to \code{\link[base]{list.files}}().}
}

\value{
  Returns an \code{\link{AffymetrixDataSet}} object.
}

\section{Detection of read maps}{
  If the first data file has a read map specified, then the corresponding
  map file is searched for and read and set to be the read map for all
  data files in this data set.  Thus, all data files must share the same
  map, if any.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
