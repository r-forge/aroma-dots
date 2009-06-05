%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSpots.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getData.BaseFileSpots}
\alias{getData.BaseFileSpots}
\alias{BaseFileSpots.getData}
\alias{getData.BaseFileSpots}
\alias{getData,BaseFileSpots-method}


 \title{Gets a subset of or the complete data table}

 \description{
   Gets a subset of or the complete data table.
 }

 \usage{\method{getData}{BaseFileSpots}(this, assays=NULL, fields=NULL, ...)}

 \arguments{
  \item{assays}{A \code{\link[base]{character}} string \code{\link[base]{vector}} of assays to be returned.}
  \item{fields}{A \code{\link[base]{character}} string \code{\link[base]{vector}} of fields, both non-assay fields ("columns") and assay fields to be returned.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{data.frame}}.
 }

 \details{
  If only one assay is returned, the assay name suffices are excluded
  from the assay field columns, otherwise they are included.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getColumns.BaseFileSpots]{*getColumns}()}
   For more information see \code{\link{BaseFileSpots}}.
 }
\keyword{internal}
\keyword{methods}
