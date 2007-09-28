%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Explorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Explorer}
\docType{class}
\alias{Explorer}


 \title{The Explorer class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class Explorer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{Explorer}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{ArrayExplorer}}, \link[aroma.affymetrix]{ChromosomeExplorer}\cr

public abstract static class \bold{Explorer}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{Explorer(tags="*", ...)}

 \arguments{
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the output path.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{addIncludes} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:display.Explorer]{display}} \tab Displays the explorer in the default browser.\cr
 \tab \code{\link[aroma.affymetrix:getAlias.Explorer]{getAlias}} \tab Gets the alias of the output set.\cr
 \tab \code{\link[aroma.affymetrix:getArrays.Explorer]{getArrays}} \tab Gets the names of the arrays.\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getName.Explorer]{getName}} \tab Gets the name of the explorer.\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.Explorer]{getTags}} \tab Gets the tags of the explorer.\cr
 \tab \code{getTagsOfInput} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.Explorer]{nbrOfArrays}} \tab Gets the total number of arrays.\cr
 \tab \code{\link[aroma.affymetrix:process.Explorer]{process}} \tab Generates image files, scripts and dynamic pages for the explorer.\cr
 \tab \code{\link[aroma.affymetrix:setAlias.Explorer]{setAlias}} \tab Sets the alias of the output set.\cr
 \tab \code{\link[aroma.affymetrix:setArrays.Explorer]{setArrays}} \tab Sets the arrays.\cr
 \tab \code{setup} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
