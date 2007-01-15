%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{importFromCrlmmFile.GenotypeCallFile}
\alias{importFromCrlmmFile.GenotypeCallFile}
\alias{GenotypeCallFile.importFromCrlmmFile}
\alias{importFromCrlmmFile.GenotypeCallFile}
\alias{importFromCrlmmFile,GenotypeCallFile-method}


\title{Imports genotype calls from a CRLMM tab-delimited call file}

\description{
 Imports genotype calls from a CRLMM tab-delimited call file.
}

\usage{\method{importFromCrlmmFile}{GenotypeCallFile}(this, filename="calls.txt", path=NULL, sample=getName(this), ..., verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the file.}
  \item{path}{An optional path to the file.}
  \item{sample}{The name of the sample column to be extracted/imported.}
  \item{...}{Not used.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns (invisibly) the number of units imported.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenotypeCallFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
