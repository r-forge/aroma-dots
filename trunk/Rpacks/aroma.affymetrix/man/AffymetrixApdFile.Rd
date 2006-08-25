%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixApdFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixApdFile}
\docType{class}
\alias{AffymetrixApdFile}


 \title{The AffymetrixApdFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixApdFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{AffymetrixApdFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixApdFile}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixDataFile}}\cr



  An AffymetrixApdFile object represents a single APD file.
 }

 \usage{AffymetrixApdFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixDataFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fromFile.AffymetrixApdFile]{fromFile}} \tab Defines an AffymetrixApdFile object from a file.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getFileType} \tab  -\cr
 \tab \code{nbrOfProbes} \tab  -\cr
 \tab \code{readIntensities} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


\bold{Methods inherited from AffymetrixDataFile}:\cr
[, addTransform, as.character, createOptimalReadMap, fromFile, getApdMap, getMapType, getProbeIntensities, getProbesetIntensities, getReadMap, hasApdMap, hasTransforms, isTransforming, normalizeQuantile, plotDensity, readCdfUnitsWriteMap, readIntensities, setApdMap, setTransform, setTransformStatus, transformAffine, transformOff, transformOn, transformProbeSignals, writeApd

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, fromFile, getFilename, getFileSize, getFileType, getName, getPath, getPathname

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
