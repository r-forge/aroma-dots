%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaMicroarrayDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFiles.AromaMicroarrayDataSet}
\alias{fromFiles.AromaMicroarrayDataSet}
\alias{AromaMicroarrayDataSet.fromFiles}
\alias{fromFiles.AromaMicroarrayDataSet}
\alias{fromFiles,AromaMicroarrayDataSet-method}


\title{Defines an AromaMicroarrayDataSet object by searching for data files}

\description{
  Defines an AromaMicroarrayDataSet object by searching for data files.
}

\synopsis{fromFiles.AromaMicroarrayDataSet(static, ..., fileClass="AromaMicroarrayDataFile", .validate=TRUE)}
\usage{AromaMicroarrayDataSet$fromFiles(..., fileClass="AromaMicroarrayDataFile", .validate=TRUE)}

\arguments{
 \item{path}{The path where to search for data files.}
 \item{...}{Optional arguments passed to the same method in the superclass.}
 \item{.validate}{If \code{\link[base:logical]{TRUE}}, the data set it validated.}
}

\value{
  Returns an \code{\link{AromaMicroarrayDataSet}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AromaMicroarrayDataSet}}.
}
\keyword{internal}
\keyword{methods}
