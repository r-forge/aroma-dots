%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readCdfBin.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{readCdfBin}
\alias{readCdfBin}


 \title{Reads a dChip CDF.bin file}

 \description{
   Reads a dChip CDF.bin file.

   Please not that this method is incomplete as it currently doesn't read
   all fields.  It is only made available so that someelse can continue
   the development.
 }

 \usage{readCdfBin(con, units=NULL, ...)}

 \arguments{
   \item{con}{A \code{\link[base:connections]{connection}} or a \code{\link[base]{character}} filename.}
   \item{units}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the units to be read.
     If \code{\link[base]{NULL}}, all units are read.}
   \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{list}} structure containing the file header and the unit data.
 }

 \examples{
path <- system.file("exData", package="dChipIO")
chipType <- "Test3"
filename <- sprintf("\%s.CDF.bin", chipType)
pathname <- file.path(path, filename)

hdr <- readCdfBinHeader(pathname)
print(hdr)

data <- readCdfBin(pathname)
str(data)

}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To read only the CDF.bin file header, see \code{\link{readCdfBinHeader}}().
 }



\keyword{file}
\keyword{IO}
