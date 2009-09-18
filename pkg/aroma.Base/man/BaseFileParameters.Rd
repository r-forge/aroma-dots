%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileParameters.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BaseFileParameters}
\docType{class}
\alias{BaseFileParameters}


 \title{The BaseFileParameters class extending BaseFileSection}

 \description{
  Package:  aroma.Base \cr
\bold{Class BaseFileParameters}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.Base]{BaseFileSection}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{BaseFileParameters}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BaseFileParameters}\cr
extends \link[aroma.Base]{BaseFileSection}\cr


 }

 \usage{BaseFileParameters(section=NULL, ...)}

 \arguments{
   \item{section}{A section \code{\link[base]{list}} structure.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.Base:as.character.BaseFileParameters]{as.character}} \tab Gets a string description of BASE parameters.\cr
 \tab \code{\link[aroma.Base:attachParameters.BaseFileParameters]{attachParameters}} \tab Assigns the parameters of a BASE parameters section locally.\cr
 \tab \code{\link[aroma.Base:getParameter.BaseFileParameters]{getParameter}} \tab Gets one parameter from a BASE parameters section.\cr
 \tab \code{\link[aroma.Base:getParameters.BaseFileParameters]{getParameters}} \tab Gets parameters from a BASE parameters section.\cr
 \tab \code{\link[aroma.Base:getPluginVersion.BaseFileParameters]{getPluginVersion}} \tab Gets the plugin version of a BASE parameters section.\cr
 \tab \code{\link[aroma.Base:nbrOfParameters.BaseFileParameters]{nbrOfParameters}} \tab Gets the number of parameters in a BASE parameters section.\cr
}


 \bold{Methods inherited from BaseFileSection}:\cr
addColumn, as.character, as.list, equals, getColumns, getData, getDataByKeys, getHeader, getHeaders, getType, hasColumns, hasData, hasHeaders, nbrOfHeaders, print, setColumns, setData, setDataFieldByKeys, setHeader, setType, sizeOfData, validate

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save

 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{BaseFile}}.
 }
\keyword{classes}
