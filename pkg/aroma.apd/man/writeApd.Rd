%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  writeApd.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{writeApd}
\alias{writeApd.default}
\alias{writeApd}

 \title{Writes an APD probe data file}

 \usage{
\method{writeApd}{default}(filename, data, ..., writeMap=NULL)
}

 \description{
   Writes an APD probe data file.
 }

 \arguments{
   \item{filename}{The filename of the APD file.}
   \item{data}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of elements to be written.}
   \item{...}{Arguments passed to \code{\link{createApd}}(), e.g. \code{chipType},
      \code{mapType} etc.}
   \item{writeMap}{A \code{\link[base]{vector}} of indicies used to change the order how
      data elements are \emph{written} (by default).}
 }

 \value{
   Returns (invisibly) the pathname to the created file.
 }

 \author{Henrik Bengtsson}

 \seealso{
   To create an APD map file, see \code{\link{readApdMap}}().
 }



\keyword{file}
\keyword{IO}
