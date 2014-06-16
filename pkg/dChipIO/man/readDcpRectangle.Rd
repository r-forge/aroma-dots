%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readDcpRectangle.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{readDcpRectangle}
\alias{readDcpRectangle}


 \title{Reads a spatial subset of probe-level data from a dChip DCP file}

 \usage{
readDcpRectangle(filename, fields=c("rawIntensities", "normalizedIntensities"),
  xrange=c(0, Inf), yrange=c(0, Inf), ..., asMatrix=TRUE)
}

 \description{
   Reads a spatial subset of probe-level data from a dChip DCP file.
 }

 \arguments{
   \item{filename}{The pathname of the DCP file.}
   \item{fields}{The cell fields to be read.}
   \item{xrange}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length two giving the left
     and right coordinates of the cells to be returned.}
   \item{yrange}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length two giving the top
     and bottom coordinates of the cells to be returned.}
   \item{asMatrix}{If \code{\link[base:logical]{TRUE}}, the CEL data fields are returned as
     matrices with element (1,1) corresponding to cell
     (xrange[1],yrange[1]).}
   \item{...}{Additional arguments passed to \code{\link{readDcp}}().}
 }

 \value{
   A named \code{\link[base]{list}} CEL structure similar to what \code{\link{readDcp}}().
   In addition, if \code{asMatrix} is \code{\link[base:logical]{TRUE}}, the CEL data fields
   are returned as matrices, otherwise not.
 }

 \author{Henrik Bengtsson}

 \examples{
path <- system.file("exData", package="dChipIO")

filename <- "Test3-1-121502.dcp"
pathname <- file.path(path, filename)

data <- readDcpRectangle(pathname)

layout(matrix(1:4, nrow=2, byrow=TRUE))
image(data$rawIntensities, main="Raw probe signals")
image(data$normalizedIntensities, main="Normalized probe signals")


}

 \seealso{
   The \code{\link{readDcp}}() method is used internally.
   This method was inspired by \code{\link[affxparser]{readCelRectangle}} of
   the \pkg{affxparser} package.
 }



\keyword{file}
\keyword{IO}
