%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFiles.GenotypeCallSet}
\alias{fromFiles.GenotypeCallSet}
\alias{GenotypeCallSet.fromFiles}
\alias{fromFiles.GenotypeCallSet}
\alias{fromFiles,GenotypeCallSet-method}


\title{Static method defining a set of genotype call files}

\description{
 Static method defining a set of genotype call files in a directory.
}

\synopsis{fromFiles.GenotypeCallSet(static, path, pattern="[.]calls$", ..., fileClass="GenotypeCallFile", verbose=FALSE)}
\usage{GenotypeCallSet$fromFiles(path, pattern="[.]calls$", ..., fileClass="GenotypeCallFile", verbose=FALSE)}

\arguments{
  \item{path}{The directory where the call files are.}
  \item{pattern}{The filename pattern used to identify call files.}
  \item{...}{Arguments passed to \code{fromFiles()} of
     \code{\link{GenericDataFileSet}}.}
  \item{fileClass}{The \code{\link[R.oo]{Class}} of each genotype file.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link{GenotypeCallSet}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenotypeCallSet}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
