%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSpots.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{setDataFields.BaseFileSpots}
\alias{setDataFields.BaseFileSpots}
\alias{BaseFileSpots.setDataFields}
\alias{setDataFields.BaseFileSpots}
\alias{setDataFields,BaseFileSpots-method}


 \title{Sets data field for an assay in a BASE file spots section}

 \description{
   Sets data field for an assay in a BASE file spots section.
 }

 \usage{\method{setDataFields}{BaseFileSpots}(this, assay=NULL, values, fields=names(values), ...)}

 \arguments{
  \item{assay}{A \code{\link[base]{character}} string for an existing assay.}
  \item{values}{A (named) \code{\link[base]{data.frame}} or a \code{\link[base]{list}} containing values to be
                                 assigned to the fields in \code{fields}.}
  \item{fields}{A \code{\link[base]{character}} of fields to be assigned.}
  \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) the new set of data as a \code{\link[base]{data.frame}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileSpots}}.
 }
\keyword{internal}
\keyword{methods}
