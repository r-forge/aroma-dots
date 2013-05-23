%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  cdfToApdMap.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{cdfToApdMap}
\alias{cdfToApdMap.default}
\alias{cdfToApdMap}

 \title{Generates an APD read map file from an Affymetrix CDF file}

 \usage{
\method{cdfToApdMap}{default}(filename, mapType=NULL, mapFile=NULL, mapPath=NULL, ..., verbose=FALSE)
}

 \description{
   Generates an APD read map file from an Affymetrix CDF file.
 }

 \arguments{
   \item{filename}{The pathname of the CDF file.}
   \item{mapType}{A \code{\link[base]{character}} string naming the type of the map.
     If \code{\link[base]{NULL}}, the chip type will be used.}
   \item{mapFile}{The filename of the resulting APD map file. If \code{\link[base]{NULL}},
     the map type with filename extension \code{apm} will be used.}
   \item{mapPath}{An optional path where to the map file will be stored.}
   \item{...}{Additional arguments passed to
      \code{\link[affxparser]{readCdfUnitsWriteMap}}, e.g. \code{units}.}
   \item{verbose}{A \code{\link[R.utils]{Verbose}} object.}
 }

 \value{
   Returns (invisibly) a \code{\link[base]{list}} structure with elements:
   \item{pathname}{The pathname of the generated APD map file.}
   \item{mapType}{The map type \code{\link[base]{character}} string.}
   \item{chipType}{The chip type \code{\link[base]{character}} string.}
   \item{readMap}{The generated read map \code{\link[base]{integer}} \code{\link[base]{vector}}.}
 }

 \author{Henrik Bengtsson}

 \seealso{
   To read an APD map file, see \code{\link{readApdMap}}().
 }



\keyword{file}
\keyword{IO}
