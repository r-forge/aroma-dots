%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{setDataFieldByKeys.BaseFileSection}
\alias{setDataFieldByKeys.BaseFileSection}
\alias{BaseFileSection.setDataFieldByKeys}
\alias{setDataFieldByKeys.BaseFileSection}
\alias{setDataFieldByKeys,BaseFileSection-method}


 \title{Sets a data field by another key field}

 \description{
   Sets a data field by another key field.
 }

 \usage{\method{setDataFieldByKeys}{BaseFileSection}(this, byField, keys, field, values, ...)}

 \arguments{
  \item{byField}{A string \code{\link[base]{character}} string of an existing column.}
  \item{keys}{Rows with these values on the key field will be returned.}
  \item{field}{A string \code{\link[base]{character}} string of an existing column.}
  \item{...}{Arguments with names of keys and values to be assigned to the
     field.}
 }

 \value{
   Returns nothing.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileSection}}.
 }
\keyword{internal}
\keyword{methods}
