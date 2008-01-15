%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.SpotData}
\alias{write.SpotData}
\alias{SpotData.write}
\alias{write.SpotData}
\alias{write,SpotData-method}


\title{Write a SpotData object to file}

\description{
  Write a SpotData object to file.
}

\usage{\method{write}{SpotData}(this, filename, path=NULL, slide=NULL, overwrite=FALSE, row.names=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the file to be written.}
  \item{path}{The path to the file.}
  \item{slide}{Index of slide to be saved.}
  \item{...}{Other arguments accepted by subclasses or which are passed
    to \code{write.table}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  spot <- SpotData$read("spot123.spot", path=system.file("data-ex", package="aroma"))

  # Write the SpotData object to a temporary file.
  filename <- paste(tempfile("SpotData"), ".dat", sep="")
  write(spot, filename)

  spot2 <- SpotData$read(filename)
  print(equals(spot, spot2))  # TRUE

  unlink(filename)
}

\seealso{
  To read one or more SpotData files at once
  see \code{\link{SpotData.read}}.
  For more information see \code{\link{SpotData}}.
}
\keyword{internal}
\keyword{methods}
