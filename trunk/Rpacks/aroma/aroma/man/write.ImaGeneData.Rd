%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ImaGeneData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.ImaGeneData}
\alias{write.ImaGeneData}
\alias{ImaGeneData.write}
\alias{write.ImaGeneData}
\alias{write,ImaGeneData-method}


\title{Write an ImaGene result data file}

\description{
 Writes the ImaGeneData object to a file using ImaGene file format.
}

\usage{\method{write}{ImaGeneData}(this, filename, path=NULL, slide=1, overwrite=FALSE, ..., digits=13, verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the file to be written.}
  \item{path}{The path to the file.}
  \item{slide}{An \code{\link[base]{integer}} specifying which slide to be written to file.
    Currently, only one slide at the time can be written.}
  \item{...}{Arguments passed to \code{write.table}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  ig <- ImaGeneData$read("imagene234.cy3", path=system.file("data-ex", package="aroma"))

  # Writes the ImaGeneData object to a file named "temp.cy3" in a format
  # that is as close as possible to the original format.
  write(ig, "temp.cy3", overwrite=TRUE)
}

\seealso{
  To read an ImaGene result data file see \code{\link[aroma:read.ImaGeneData]{*read}()}.
  For more information see \code{\link{ImaGeneData}}.
}
\keyword{internal}
\keyword{methods}
