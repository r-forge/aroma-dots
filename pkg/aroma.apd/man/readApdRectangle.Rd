%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readApdRectangle.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readApdRectangle}
\alias{readApdRectangle.default}
\alias{readApdRectangle}

\title{Reads a spatial subset of probe-level data from Affymetrix APD files}

\usage{
\method{readApdRectangle}{default}(filename, xrange=c(0, Inf), yrange=c(0, Inf), ..., asMatrix=TRUE)
}

\description{
  Reads a spatial subset of probe-level data from Affymetrix APD files.
}

\arguments{
  \item{filename}{The pathname of the APD file.}
  \item{xrange}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length two giving the left
    and right coordinates of the cells to be returned.}
  \item{yrange}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length two giving the top
    and bottom coordinates of the cells to be returned.}
  \item{...}{Additional arguments passed to \code{\link{readApd}}().}
  \item{asMatrix}{If \code{\link[base:logical]{TRUE}}, the APD data fields are returned as
    matrices with element (1,1) corresponding to cell
    (xrange[1],yrange[1]).}
}

\value{
  A named \code{\link[base]{list}} APD structure similar to what \code{\link{readApd}}().
  In addition, if \code{asMatrix} is \code{\link[base:logical]{TRUE}}, the APD data fields
  are returned as matrices, otherwise not.
}

\author{Henrik Bengtsson}

\examples{

# Local functions
rotate270 <- function(x, ...) {
  x <- t(x)
  nc <- ncol(x)
  if (nc < 2) return(x)
  x[,nc:1,drop=FALSE]
}

# Scan current directory for APD files
files <- list.files(pattern="[.](apd|APD)$")
files <- files[!file.info(files)$isdir]
if (length(files) > 0) {
  apdFile <- files[1]
 
  # Read APD intensities in the upper left corner
  apd <- readApdRectangle(apdFile, xrange=c(0,250), yrange=c(0,250))
  z <- rotate270(apd$intensities)
  sub <- paste("Chip type:", apd$header$chipType)
  image(z, col=gray.colors(256), axes=FALSE, main=apdFile, sub=sub)
  text(x=0, y=1, labels="(0,0)", adj=c(0,-0.7), cex=0.8, xpd=TRUE)
  text(x=1, y=0, labels="(250,250)", adj=c(1,1.2), cex=0.8, xpd=TRUE)
}
}

\seealso{
  The \code{\link{readApd}}() method is used internally.
}



\keyword{file}
\keyword{IO}
