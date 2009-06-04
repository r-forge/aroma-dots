%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readDcp.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{readDcp}
\alias{readDcp}


 \title{Reads a dChip DCP file}

 \description{
   Reads a dChip DCP file.
 }

 \usage{readDcp(con, fields=c("rawIntensities", "normalizedIntensities", "calls", "thetas", "thetaStds", "excludes"), cells=NULL, units=NULL, .nbrOfUnits=NULL, ...)}

 \arguments{
   \item{con}{A \code{\link[base:connections]{connection}} or a \code{\link[base]{character}} filename.}
   \item{fields}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying the fields to be read.}
   \item{cells}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the indices of the cell
     data to be read.}
   \item{units}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the indices of the unit
     data to be read.}
   \item{.nbrOfUnits}{A \code{\link[base]{integer}} specifying the number of units available
     in the file. If \code{\link[base]{NULL}}, this is inferred from the file size and the
     file header. The dChip software itself instead uses the corrsponding
     value in the CDF.bin file, but that file is specified by the user
     leaving room for errors.}
   \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{list}} structure containing the file header and the
   requested data fields.
 }

 \examples{
path <- system.file("exData", package="dChipIO")

filename <- "Test3-1-121502.dcp"
pathname <- file.path(path, filename)

hdr <- readDcpHeader(pathname)
print(hdr)

data <- readDcp(pathname)
str(data)

}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To read only the DCP file header, see \code{\link{readDcpHeader}}().
 }



\keyword{file}
\keyword{IO}
