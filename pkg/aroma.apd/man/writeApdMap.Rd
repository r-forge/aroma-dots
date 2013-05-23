%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  writeApdMap.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{writeApdMap}
\alias{writeApdMap.default}
\alias{writeApdMap}

 \title{Writes an APD probe map file}

 \usage{
\method{writeApdMap}{default}(filename, path=NULL, map, ...)
}

 \description{
   Writes an APD probe map file.
 }

 \arguments{
   \item{filename}{The filename of the APD file.}
   \item{path}{The path to the APD file.}
   \item{map}{A \code{\link[base]{vector}} of indicies.}
   \item{...}{Additional arguments passed to \code{\link{writeApd}}().}
 }

 \value{
   Returns (invisibly) the pathname to the create file.
 }

 \author{Henrik Bengtsson}

 \seealso{
   To read an APD map file, see \code{\link{readApdMap}}().
 }



\keyword{file}
\keyword{IO}
