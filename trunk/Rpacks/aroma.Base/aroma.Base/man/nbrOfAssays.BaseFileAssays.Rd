%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileAssays.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{nbrOfAssays.BaseFileAssays}
\alias{nbrOfAssays.BaseFileAssays}
\alias{BaseFileAssays.nbrOfAssays}
\alias{nbrOfAssays.BaseFileAssays}
\alias{nbrOfAssays,BaseFileAssays-method}


 \title{Gets the (expected) number of assays in the data table}

 \description{
   Gets the (expected) number of assays in the data table.  This method return what \code{\link[aroma.Base:getCount.BaseFileAssays]{*getCount}()} returns,
   except that if the latter returns \code{\link[base]{NULL}}, this returns zero.
 }

 \usage{\method{nbrOfAssays}{BaseFileAssays}(this, ...)}

 \arguments{
  \item{...}{Not used.}
 }

 \value{
   Returns an \code{\link[base]{integer}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getCount.BaseFileAssays]{*getCount}()}
   For more information see \code{\link{BaseFileAssays}}.
 }
\keyword{internal}
\keyword{methods}
