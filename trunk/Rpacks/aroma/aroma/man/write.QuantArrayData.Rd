%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantArrayData.R
% 
% on Tue Jan 15 18:35:00 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.QuantArrayData}
\alias{write.QuantArrayData}
\alias{QuantArrayData.write}
\alias{write.QuantArrayData}
\alias{write,QuantArrayData-method}


\title{Write a QuantArray result data file}

\description{
 Writes the QuantArrayData object to a file using QuantArray file format.
}

\usage{\method{write}{QuantArrayData}(this, filename, path=NULL, slide=1, overwrite=FALSE, ..., digits=9, verbose=FALSE)}

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
 \dontrun{
  # At the moment there is no QuantArray sample files in the package...
  # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
  # ...will use a GenePix sample file instead.
  qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  # Writes the QuantArrayData object to a file named "temp.txt" in a format
  # that is as close as possible to the original format.
  write(qa, "temp.txt", overwrite=TRUE)
 }
}

\seealso{
  To read a QuantArray result data file see \code{\link[aroma:read.QuantArrayData]{*read}()}.
  For more information see \code{\link{QuantArrayData}}.
}
\keyword{internal}
\keyword{methods}
