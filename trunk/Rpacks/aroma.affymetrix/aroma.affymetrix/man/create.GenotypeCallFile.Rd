%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{create.GenotypeCallFile}
\alias{create.GenotypeCallFile}
\alias{GenotypeCallFile.create}
\alias{create.GenotypeCallFile}
\alias{create,GenotypeCallFile-method}


\title{Static method creating a genotype call file}

\description{
 Static method creating a genotype call file.
}

\synopsis{create.GenotypeCallFile(static, filename, path=NULL, cdf, overwrite=FALSE, ...)}
\usage{GenotypeCallFile$create(filename, path=NULL, cdf, overwrite=FALSE, ...)}

\arguments{
  \item{filename}{The filename of the file.}
  \item{path}{An optional path to the file.}
  \item{cdf}{An \code{\link{AffymetrixCdfFile}}.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, an already existing file is overwritten,
     otherwise not.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link{GenotypeCallFile}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenotypeCallFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
