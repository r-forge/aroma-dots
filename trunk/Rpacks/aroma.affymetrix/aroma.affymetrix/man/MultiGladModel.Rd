%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiGladModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MultiGladModel}
\docType{class}
\alias{MultiGladModel}


 \title{The MultiGladModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class MultiGladModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{MultiGladModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{MultiGladModel}\cr
extends \link[R.oo]{Object}\cr



  This class is similar to the \code{\link{GladModel}} class, but it can handle
  multiple chip types.  For now, consider this class to be very much
  a trial-and-error class.
 }

 \usage{MultiGladModel(cesList=NULL, referenceList=NULL, tags="*", ...)}

 \arguments{
   \item{cesList}{A \code{\link[base]{list}} of \code{\link{ChipEffectSet}} objects.}
   \item{referenceList}{A \code{\link{ChipEffectFile}} objects.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{fit} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fitOne.MultiGladModel]{fitOne}} \tab Fits the GLAD model for one chromosome in one sample.\cr
 \tab \code{getArrays} \tab  -\cr
 \tab \code{getChipEffectFiles} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.MultiGladModel]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getListOfChipEffects} \tab  -\cr
 \tab \code{getListOfGenomeInformations} \tab  -\cr
 \tab \code{getListOfReferences} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getReferenceFiles} \tab  -\cr
 \tab \code{getReferenceName} \tab  -\cr
 \tab \code{getRegions} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTableOfArrays} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{nbrOfChipTypes} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{writeRegions} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  \code{\link{GladModel}}.
 }
\keyword{classes}
