%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSpots.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BaseFileSpots}
\docType{class}
\alias{BaseFileSpots}


 \title{The BaseFileSpots class extending BaseFileSection}

 \description{
  Package:  aroma.Base \cr
\bold{Class BaseFileSpots}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.Base]{BaseFileSection}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{BaseFileSpots}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BaseFileSpots}\cr
extends \link[aroma.Base]{BaseFileSection}\cr


 }

 \usage{BaseFileSpots(...)}

 \arguments{
   \item{...}{Arguments passed to constructor of superclass.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.Base:as.character.BaseFileSpots]{as.character}} \tab Gets a string description of object.\cr
 \tab \code{\link[aroma.Base:getAssayFields.BaseFileSpots]{getAssayFields}} \tab Gets the names of the assay-field data columns.\cr
 \tab \code{\link[aroma.Base:getAssays.BaseFileSpots]{getAssays}} \tab Gets the ids of the assays in this section.\cr
 \tab \code{\link[aroma.Base:getColumns.BaseFileSpots]{getColumns}} \tab Gets the names of the non-assay-field data columns.\cr
 \tab \code{\link[aroma.Base:getData.BaseFileSpots]{getData}} \tab Gets a subset of or the complete data table.\cr
 \tab \code{\link[aroma.Base:getDataFiles.BaseFileSpots]{getDataFiles}} \tab Gets the filenames where data of an assay-field section is stored.\cr
 \tab \code{\link[aroma.Base:hasDataFiles.BaseFileSpots]{hasDataFiles}} \tab Checks if data of an assay-field section is stored on file.\cr
 \tab \code{\link[aroma.Base:nbrOfAssays.BaseFileSpots]{nbrOfAssays}} \tab Gets the number of assays in this section.\cr
 \tab \code{\link[aroma.Base:setAssays.BaseFileSpots]{setAssays}} \tab Sets the ids of the assays in this section.\cr
 \tab \code{\link[aroma.Base:setColumns.BaseFileSpots]{setColumns}} \tab Sets the column names of the data table.\cr
 \tab \code{\link[aroma.Base:setDataFields.BaseFileSpots]{setDataFields}} \tab Sets data field for an assay in a BASE file spots section.\cr
 \tab \code{\link[aroma.Base:setDataFiles.BaseFileSpots]{setDataFiles}} \tab Sets the filenames where data of an assay-field section is stored.\cr
}


 \bold{Methods inherited from BaseFileSection}:\cr
addColumn, as.character, as.list, equals, getColumns, getData, getDataByKeys, getHeader, getHeaders, getType, hasColumns, hasData, hasHeaders, nbrOfHeaders, print, setColumns, setData, setDataFieldByKeys, setHeader, setType, sizeOfData, validate

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
