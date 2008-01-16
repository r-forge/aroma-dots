%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotfinderData.R
% 
% on Tue Jan 15 18:37:29 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.SpotfinderData}
\alias{write.SpotfinderData}
\alias{SpotfinderData.write}
\alias{write.SpotfinderData}
\alias{write,SpotfinderData-method}


\title{Write a SpotfinderData object to file}

\description{
 Write a SpotfinderData object to file.
}

\usage{\method{write}{SpotfinderData}(this, filename, path=NULL, slide=NULL, overwrite=FALSE, verbose=FALSE, ...)}

\arguments{
  \item{filename}{The filename of the file to be written.}
  \item{path}{The path to the file.}
  \item{slide}{Slide to be saved.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{
  sf <- SpotfinderData$read("sf111.tav", path=system.file("data-ex", package="aroma"))

  # Write the SpotfinderData object to a temporary file.
  filename <- paste(tempfile("SpotfinderData"), ".tav", sep="")
  write(sf, filename)

  sf2 <- SpotfinderData$read(filename)
  print(equals(sf, sf2))  # TRUE

  unlink(filename)
}}

\seealso{
  To read one or more SpotfinderData files at once
  see \code{\link[aroma:read.SpotfinderData]{*read}()}.
  For more information see \code{\link{SpotfinderData}}.
}
\keyword{internal}
\keyword{methods}
