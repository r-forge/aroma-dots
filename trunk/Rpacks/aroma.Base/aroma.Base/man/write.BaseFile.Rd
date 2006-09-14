%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{write.BaseFile}
\alias{write.BaseFile}
\alias{BaseFile.write}
\alias{write.BaseFile}
\alias{write,BaseFile-method}


 \title{Writes this BASE file structure}

 \description{
   Writes this BASE file structure.
 }

 \usage{\method{write}{BaseFile}(this, file, ...)}

 \arguments{
  \item{file}{A filename or a \code{\link[base:connections]{connection}} to be written.}
  \item{...}{If \code{file} is not a \code{\link[base:connections]{connection}}, arguments passed to
    to \code{\link[R.utils]{getWritablePathname.Arguments}()},
    e.g. \code{path}, \code{overwrite}, and \code{mkdirs}.}
 }

 \value{
   Returns what \code{\link{writeBaseFile}}() returns.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:read.BaseFile]{*read}()}.
   For more information see \code{\link{BaseFile}}.
 }
\keyword{internal}
\keyword{methods}
