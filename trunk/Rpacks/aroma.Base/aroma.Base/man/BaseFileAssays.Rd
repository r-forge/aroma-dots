%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileAssays.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BaseFileAssays}
\docType{class}
\alias{BaseFileAssays}


 \title{The BaseFileAssays class extending BaseFileSection}

 \description{
  Package:  aroma.Base \cr
\bold{Class BaseFileAssays}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.Base]{BaseFileSection}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{BaseFileAssays}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BaseFileAssays}\cr
extends \link[aroma.Base]{BaseFileSection}\cr


 }

 \usage{BaseFileAssays(...)}

 \arguments{
   \item{...}{Arguments passed to constructor of superclass.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.Base:as.character.BaseFileAssays]{as.character}} \tab Gets a string description of object.\cr
 \tab \code{\link[aroma.Base:getAnnotationColumns.BaseFileAssays]{getAnnotationColumns}} \tab Gets the names of the annotation columns.\cr
 \tab \code{\link[aroma.Base:getCount.BaseFileAssays]{getCount}} \tab Gets the (expected) number of assays in the data table.\cr
 \tab \code{\link[aroma.Base:getDataByIds.BaseFileAssays]{getDataByIds}} \tab Gets a subset of data for one or all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:getId.BaseFileAssays]{getId}} \tab Gets the id of all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:getName.BaseFileAssays]{getName}} \tab Gets the name of all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:getParents.BaseFileAssays]{getParents}} \tab Gets the ids of the assay parents of all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:nbrOfAssays.BaseFileAssays]{nbrOfAssays}} \tab Gets the (expected) number of assays in the data table.\cr
 \tab \code{\link[aroma.Base:setAnnotationColumns.BaseFileAssays]{setAnnotationColumns}} \tab Sets the names of the annotation columns.\cr
 \tab \code{\link[aroma.Base:setCount.BaseFileAssays]{setCount}} \tab Sets (or update) the number of assays in the data table.\cr
 \tab \code{\link[aroma.Base:setData.BaseFileAssays]{setData}} \tab Sets data of a BASE file section.\cr
 \tab \code{\link[aroma.Base:setDataFieldByIds.BaseFileAssays]{setDataFieldByIds}} \tab Sets a field of one or all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:setName.BaseFileAssays]{setName}} \tab Sets the names of one or all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:setParents.BaseFileAssays]{setParents}} \tab Sets the ids of the assay parents of one or all bioassays in the bioassay set.\cr
 \tab \code{\link[aroma.Base:validate.BaseFileAssays]{validate}} \tab Validates assays section.\cr
}


\bold{Methods inherited from BaseFileSection}:\cr
addColumn, as.character, as.list, equals, getColumns, getData, getDataByKeys, getHeader, getHeaders, getType, hasColumns, hasData, hasHeaders, nbrOfHeaders, print, setColumns, setData, setDataFieldByKeys, setHeader, setType, sizeOfData, validate

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, equals, extend, finalize, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
