%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BaseFileSection}
\docType{class}
\alias{BaseFileSection}


 \title{The BaseFileSection class}

 \description{
  Package:  aroma.Base \cr
\bold{Class BaseFileSection}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{BaseFileSection}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.Base]{BaseFileAssays}, \link[aroma.Base]{BaseFileBioassaysetinfo}, \link[aroma.Base]{BaseFileParameters}, \link[aroma.Base]{BaseFileSpots}\cr

public static class \bold{BaseFileSection}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{BaseFileSection(section=NULL, type=NULL, ...)}

 \arguments{
   \item{section}{A \code{\link[base]{list}} structure.}
   \item{type}{A \code{\link[base]{character}} string specifying the type of section.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.Base:addColumn.BaseFileSection]{addColumn}} \tab Adds a column to the other column names.\cr
 \tab \code{\link[aroma.Base:as.character.BaseFileSection]{as.character}} \tab Gets a string description of object.\cr
 \tab \code{\link[aroma.Base:as.list.BaseFileSection]{as.list}} \tab Gets a list representation of object.\cr
 \tab \code{\link[aroma.Base:equals.BaseFileSection]{equals}} \tab Checks if a BASE file section equals another.\cr
 \tab \code{\link[aroma.Base:getColumns.BaseFileSection]{getColumns}} \tab Gets the column names of the data table.\cr
 \tab \code{\link[aroma.Base:getData.BaseFileSection]{getData}} \tab Gets all or a subset of data of a BASE file section.\cr
 \tab \code{\link[aroma.Base:getDataByKeys.BaseFileSection]{getDataByKeys}} \tab Gets a subset of data for one or all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:getHeader.BaseFileSection]{getHeader}} \tab Gets one header of a BASE section.\cr
 \tab \code{\link[aroma.Base:getHeaders.BaseFileSection]{getHeaders}} \tab Gets the headers of a BASE section.\cr
 \tab \code{\link[aroma.Base:getType.BaseFileSection]{getType}} \tab Gets the type of a BASE section.\cr
 \tab \code{\link[aroma.Base:hasColumns.BaseFileSection]{hasColumns}} \tab Checks if certain columns exists or not.\cr
 \tab \code{\link[aroma.Base:hasData.BaseFileSection]{hasData}} \tab Checks if section has a data table.\cr
 \tab \code{\link[aroma.Base:hasHeaders.BaseFileSection]{hasHeaders}} \tab Checks if specified headers exists in a BASE section.\cr
 \tab \code{\link[aroma.Base:nbrOfHeaders.BaseFileSection]{nbrOfHeaders}} \tab Gets the number of headers in a BASE section.\cr
 \tab \code{\link[aroma.Base:print.BaseFileSection]{print}} \tab Prints a BASE file section to stdout.\cr
 \tab \code{\link[aroma.Base:setColumns.BaseFileSection]{setColumns}} \tab Sets the column names of the data table.\cr
 \tab \code{\link[aroma.Base:setData.BaseFileSection]{setData}} \tab Sets data of a BASE file section.\cr
 \tab \code{\link[aroma.Base:setDataFieldByKeys.BaseFileSection]{setDataFieldByKeys}} \tab Sets a data field by another key field.\cr
 \tab \code{\link[aroma.Base:setHeader.BaseFileSection]{setHeader}} \tab Sets a header of a BASE section.\cr
 \tab \code{\link[aroma.Base:setType.BaseFileSection]{setType}} \tab Sets the type of a BASE section.\cr
 \tab \code{\link[aroma.Base:sizeOfData.BaseFileSection]{sizeOfData}} \tab Gets the dimension of the data table.\cr
 \tab \code{\link[aroma.Base:validate.BaseFileSection]{validate}} \tab Validates BASE file section.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, equals, extend, finalize, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

 }

 \details{
 }


 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{BaseFile}}.
 }
\keyword{classes}
