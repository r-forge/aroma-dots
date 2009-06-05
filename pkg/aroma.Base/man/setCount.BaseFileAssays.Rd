%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileAssays.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{setCount.BaseFileAssays}
\alias{setCount.BaseFileAssays}
\alias{BaseFileAssays.setCount}
\alias{setCount.BaseFileAssays}
\alias{setCount,BaseFileAssays-method}


 \title{Sets (or update) the number of assays in the data table}

 \description{
   Sets (or update) the number of assays in the data table.
 }

 \usage{\method{setCount}{BaseFileAssays}(this, count=NULL, ...)}

 \arguments{
  \item{count}{An \code{\link[base]{integer}} or \code{\link[base]{NULL}}. If \code{\link[base]{NULL}}, the number of assays is
    set to the number of rows in the current data table.}
  \item{...}{Not used.}
 }

 \value{
   Returns an \code{\link[base]{integer}} for the old count, or \code{\link[base]{NULL}} if missing.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getCount.BaseFileAssays]{*getCount}()}
   For more information see \code{\link{BaseFileAssays}}.
 }
\keyword{internal}
\keyword{methods}
