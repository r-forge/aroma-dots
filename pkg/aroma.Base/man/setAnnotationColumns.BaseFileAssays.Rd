%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileAssays.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{setAnnotationColumns.BaseFileAssays}
\alias{setAnnotationColumns.BaseFileAssays}
\alias{BaseFileAssays.setAnnotationColumns}
\alias{setAnnotationColumns.BaseFileAssays}
\alias{setAnnotationColumns,BaseFileAssays-method}


 \title{Sets the names of the annotation columns}

 \description{
   Sets the names of the annotation columns.
 }

 \usage{\method{setAnnotationColumns}{BaseFileAssays}(this, names, ...)}

 \arguments{
  \item{names}{A \code{\link[base]{character}} \code{\link[base]{vector}} of column names.
    If there already are existing annotation columns specified, the
    number of new column names must match the existing number, otherwise
    an exception is thrown.}
  \item{...}{Not used.}
 }

 \value{
   Returns nothing.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getAnnotationColumns.BaseFileAssays]{*getAnnotationColumns}()}
   For more information see \code{\link{BaseFileAssays}}.
 }
\keyword{internal}
\keyword{methods}
