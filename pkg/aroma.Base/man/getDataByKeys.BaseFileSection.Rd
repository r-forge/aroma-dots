%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getDataByKeys.BaseFileSection}
\alias{getDataByKeys.BaseFileSection}
\alias{BaseFileSection.getDataByKeys}
\alias{getDataByKeys.BaseFileSection}
\alias{getDataByKeys,BaseFileSection-method}


 \title{Gets a subset of data for one or all bioassays in the bioassay set}

 \description{
   Gets a subset of data for one or all bioassays in the bioassay set.
 }

 \usage{\method{getDataByKeys}{BaseFileSection}(this, byField=NULL, keys=NULL, fields=NULL, ...)}

 \arguments{
  \item{byField}{A \code{\link[base]{character}} string of the field to be stratified on.
     If \code{\link[base]{NULL}}, no stratification is done.}
  \item{keys}{Rows with these values on the key field will be returned.}
  \item{fields}{A \code{\link[base]{character}} \code{\link[base]{vector}} of the fields to be returned.
     If \code{\link[base]{NULL}}, all fields are returned.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{data.frame}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getData.BaseFileSection]{*getData}()}.
   For more information see \code{\link{BaseFileSection}}.
 }
\keyword{internal}
\keyword{methods}
