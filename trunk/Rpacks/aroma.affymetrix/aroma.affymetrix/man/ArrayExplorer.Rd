%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ArrayExplorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ArrayExplorer}
\docType{class}
\alias{ArrayExplorer}


 \title{The ArrayExplorer class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ArrayExplorer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{ArrayExplorer}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ArrayExplorer}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{ArrayExplorer(celSet=NULL, tags="*", ...)}

 \arguments{
   \item{celSet}{An \code{\link{AffymetrixCelSet}} object.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the output path.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{addColorMap} \tab  -\cr
 \tab \code{addIncludes} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:display.ArrayExplorer]{display}} \tab Displays the explorer in the default browser.\cr
 \tab \code{\link[aroma.affymetrix:getAlias.ArrayExplorer]{getAlias}} \tab Gets the alias of the output set.\cr
 \tab \code{\link[aroma.affymetrix:getArrays.ArrayExplorer]{getArrays}} \tab Gets the names of the arrays.\cr
 \tab \code{getColorMaps} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getDataSet.ArrayExplorer]{getDataSet}} \tab Gets the data set.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getName.ArrayExplorer]{getName}} \tab Gets the name of the explorer.\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.ArrayExplorer]{getTags}} \tab Gets the tags of the explorer.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.ArrayExplorer]{nbrOfArrays}} \tab Gets the total number of arrays.\cr
 \tab \code{\link[aroma.affymetrix:process.ArrayExplorer]{process}} \tab Generates image files, scripts and dynamic pages for the explorer.\cr
 \tab \code{\link[aroma.affymetrix:setAlias.ArrayExplorer]{setAlias}} \tab Sets the alias of the output set.\cr
 \tab \code{\link[aroma.affymetrix:setArrays.ArrayExplorer]{setArrays}} \tab Sets the arrays.\cr
 \tab \code{setColorMaps} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
