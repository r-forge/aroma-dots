%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSetTuple.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixCelSetTuple}
\docType{class}
\alias{AffymetrixCelSetTuple}


 \title{The AffymetrixCelSetTuple class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCelSetTuple}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{AffymetrixCelSetTuple}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ChipEffectSetTuple}\cr

public static class \bold{AffymetrixCelSetTuple}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{AffymetrixCelSetTuple(csList=NULL, tags="*", ..., .setClass="AffymetrixCelSet")}

 \arguments{
   \item{csList}{A single or \code{\link[base]{list}} of \code{\link{AffymetrixCelSet}}:s.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{...}{Not used.}
   \item{.setClass}{The name of the class of the input set.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byPath} \tab  -\cr
 \tab \code{clearCache} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:extract.AffymetrixCelSetTuple]{extract}} \tab Extracts a subset AffymetrixCelSetTuple.\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getArrays.AffymetrixCelSetTuple]{getArrays}} \tab Gets the names of the arrays.\cr
 \tab \code{\link[aroma.affymetrix:getArrayTuple.AffymetrixCelSetTuple]{getArrayTuple}} \tab Gets arrays across chip types for one sample.\cr
 \tab \code{\link[aroma.affymetrix:getChipTypes.AffymetrixCelSetTuple]{getChipTypes}} \tab Gets a label for all chip types merged.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{getListOfSets} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.AffymetrixCelSetTuple]{nbrOfArrays}} \tab Gets the number of arrays.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfChipTypes.AffymetrixCelSetTuple]{nbrOfChipTypes}} \tab Gets the number of chip types.\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
