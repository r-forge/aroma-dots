%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readApdMap.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{readApdMap}
\alias{readApdMap.default}
\alias{readApdMap}

 \title{Reads an APD probe map file}

 \usage{\method{readApdMap}{default}(filename, path=NULL, ...)}

 \description{
   Reads an APD probe map file.
 }

 \arguments{
   \item{filename}{The filename of the APD file.}
   \item{path}{The path to the APD file.}
   \item{...}{Arguments passed to \code{\link{readApd}}().}
 }

 \value{
   A named \code{\link[base]{list}} with the two elements \code{header} and
   \code{map}.  The header is in turn a \code{\link[base]{list}} structure and
   the second is a \code{\link[base]{numeric}} \code{\link[base]{vector}} holding the probe map indices.
 }

 \section{File format}{
   The file format of an APD map file is identical to the file format
   of an APD file, see \code{\link{readApd}}().  The APD map file identified by
   the name of the data defaults to \code{"map"}.  If not, an error
   is thrown.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To search for an APD map file, see \code{\link{findApdMap}}().
   To create a cell index map from an CDF file, see
   \code{\link[affxparser]{readCdfUnitsMap}}.
   Internally, \code{\link{readApd}}() is used.
 }



\keyword{file}
\keyword{IO}
