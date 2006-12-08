%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PlatformDesign.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PlatformDesign}
\docType{class}
\alias{PlatformDesign}


 \title{The PlatformDesign class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class PlatformDesign}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{PlatformDesign}\cr

\bold{Directly known subclasses:}\cr
\cr

public class \bold{PlatformDesign}\cr
extends \link[R.oo]{Object}\cr



  The PlatformDesign class provides a memory-efficient interface to a
  so called platform-design package.  Its method can access a subset of the
  data of such packages without having to load the package, which saves
  memory.
 }

 \usage{PlatformDesign(cdf=NULL, ...)}

 \arguments{
  \item{cdf}{An \code{\link{AffymetrixCdfFile}}.}
  \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{clearFileCache} \tab  -\cr
 \tab \code{getAnnotations} \tab  -\cr
 \tab \code{getCachedFields} \tab  -\cr
 \tab \code{getCachePath} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getCrlmmInfo} \tab  -\cr
 \tab \code{getDataFile} \tab  -\cr
 \tab \code{getDataFiles} \tab  -\cr
 \tab \code{getDescription} \tab  -\cr
 \tab \code{getFeatureInfo} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPackageName} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getPlatformDesignObject} \tab  -\cr
 \tab \code{getProbeEffects} \tab  -\cr
 \tab \code{getReferenceQuantiles} \tab  -\cr
 \tab \code{isInstalled} \tab  -\cr
 \tab \code{isLoaded} \tab  -\cr
 \tab \code{load} \tab  -\cr
 \tab \code{unload} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
