%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ApdMap.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ApdMap}
\docType{class}
\alias{ApdMap}


 \title{The ApdMap class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ApdMap}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{ApdMap}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ApdMap}\cr
extends \link[R.oo]{Object}\cr



  An ApdMap object represents an APD read and write map.
  Probe indices in an APD data file may be reorder, for instance, such
  that elements are read in contiguous blocks instead of in an randomized
  order which is the typical order in CEL files.  This will speed up reading.
  An ApdMap object keeps track of the index maps for reading and writing
  such files.
 }

 \usage{ApdMap(chipType=NULL, mapType=chipType, ...)}

 \arguments{
   \item{chipType}{A \code{\link[base]{character}} string specifying the chip type.}
   \item{mapType}{A \code{\link[base]{character}} string specifying the map type.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{equals} \tab  -\cr
 \tab \code{findApdMap} \tab  -\cr
 \tab \code{fromCdf} \tab  -\cr
 \tab \code{fromMapType} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getMapType} \tab  -\cr
 \tab \code{getReadMap} \tab  -\cr
 \tab \code{getWriteMap} \tab  -\cr
 \tab \code{read} \tab  -\cr
 \tab \code{write} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \seealso{
  \code{\link[aroma.apd]{readApdMap}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
