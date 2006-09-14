%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{equals.BaseFile}
\alias{equals.BaseFile}
\alias{BaseFile.equals}
\alias{equals.BaseFile}
\alias{equals,BaseFile-method}


 \title{Checks if this object equals another}

 \description{
   Checks if this object equals another.
 }

 \usage{\method{equals}{BaseFile}(this, other, ...)}

 \arguments{
  \item{other}{Another object.}
  \item{...}{Not used.}
 }

 \value{
   Returns \code{\link[base:logical]{TRUE}} if the objects are equal, otherwise \code{\link[base:logical]{FALSE}}.
 }

 \details{
  If the other object is not a BaseFile object, this method returns \code{\link[base:logical]{FALSE}}.
  Otherwise, \code{as.list()} is called on both objects are these are
  compared with \code{\link[base]{identical}}().
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFile}}.
 }
\keyword{internal}
\keyword{methods}
