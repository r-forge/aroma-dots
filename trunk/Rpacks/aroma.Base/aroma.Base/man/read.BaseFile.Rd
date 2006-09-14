%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{read.BaseFile}
\alias{read.BaseFile}
\alias{BaseFile.read}
\alias{read.BaseFile}
\alias{read,BaseFile-method}


 \title{Static method to read a BASE file structure}

 \description{
   Static method to read a BASE file structure.
 }

 \synopsis{read.BaseFile(static, file, path=NULL, ...)}
\usage{BaseFile$read(file, path=NULL, ...)}

 \arguments{
  \item{file}{A filename or a \code{\link[base:connections]{connection}} to be read.}
  \item{path}{If \code{file} is not a \code{\link[base:connections]{connection}}, an optional path to
    be added to the file to be read.}
  \item{...}{Additional arguments passed to \code{\link{readBaseFile}}(),
    especially \code{extractSpotsData}.}
 }

 \value{
   Returns an \code{\link{BaseFile}} object.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:write.BaseFile]{*write}()}.
   For more information see \code{\link{BaseFile}}.
 }
\keyword{internal}
\keyword{methods}
