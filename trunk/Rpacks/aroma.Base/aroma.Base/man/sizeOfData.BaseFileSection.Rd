%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{sizeOfData.BaseFileSection}
\alias{sizeOfData.BaseFileSection}
\alias{BaseFileSection.sizeOfData}
\alias{sizeOfData.BaseFileSection}
\alias{sizeOfData,BaseFileSection-method}


 \title{Gets the dimension of the data table}

 \description{
   Gets the dimension of the data table, if it exists.
 }

 \usage{\method{sizeOfData}{BaseFileSection}(this, ...)}

 \arguments{
  \item{...}{Not used.}
 }

 \value{
   If data table exists, an \code{\link[base]{integer}} \code{\link[base]{vector}} of length two specifying the
   number of rows and columns of the table is returned. Otherwise, \code{\link[base]{NULL}}
   is returned.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileSection}}.
 }
\keyword{internal}
\keyword{methods}
