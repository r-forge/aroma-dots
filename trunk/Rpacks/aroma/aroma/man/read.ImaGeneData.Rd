%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ImaGeneData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.ImaGeneData}
\alias{read.ImaGeneData}
\alias{ImaGeneData.read}
\alias{read.ImaGeneData}
\alias{read,ImaGeneData-method}


\title{Reads an ImaGene result file}

\description{
 Reads an ImaGene result file.
 Currently ImaGene v4.1, v4.2, v5.0 files are supported. If the version
 is not recognized it will try to make a best guess and read it in anyway.
 However, it can not be guaranteed that all fields will be of the correct
 data type.
}

\synopsis{read.ImaGeneData(this, filename, path=NULL, verbose=FALSE, ...)}
\usage{ImaGeneData$read(filename, path=NULL, verbose=FALSE, ...)}

\arguments{
  \item{filename}{The filename of the file to be read.}
  \item{path}{The path to the file.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns an \code{\link{ImaGeneData}} object.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  igG <- ImaGeneData$read("imagene234.cy3", path=system.file("data-ex", package="aroma"))
  igR <- ImaGeneData$read("imagene234.cy5", path=system.file("data-ex", package="aroma"))
  raw <- getRawData(igR, igG)
}

\seealso{
  For more information see \code{\link{ImaGeneData}}.
}
\keyword{internal}
\keyword{methods}
