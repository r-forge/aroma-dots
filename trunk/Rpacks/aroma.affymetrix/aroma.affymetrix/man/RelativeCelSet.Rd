%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RelativeCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RelativeCelSet}
\docType{class}
\alias{RelativeCelSet}


 \title{The RelativeCelSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class RelativeCelSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{RelativeCelSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{RelativeCelSet}\cr
extends \link[aroma.affymetrix]{AffymetrixCelSet}\cr



  This class represents an AffymetrixCelSet whose signals are calculated
  relative to the signals of a reference AffymetrixCelFile.
 }

 \usage{RelativeCelSet(dataSet=NULL, reference=NULL, relativeFields=c("intensities"), ...)}

 \arguments{
   \item{dataSet}{A \code{\link{AffymetrixCelSet}}.}
  \item{...}{Arguments passed to the constructor of \code{\link{AffymetrixCelSet}}.}
  \item{reference}{A reference \code{\link{AffymetrixCelFile}} object with a
    compatible chip type.}
  \item{relativeFields}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying for which fields
    relative signals should be calculated.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getData} \tab  -\cr
 \tab \code{getReferenceFile} \tab  -\cr
 \tab \code{getRelativeFields} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


\bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, append, applyToUnitIntensities, as, as.AffymetrixCelSet, as.character, averageQuantile, bgAdjustGcrma, bgAdjustOptical, bgAdjustRma, calculateGsbParameters, calibrateAllelicCrosstalk, clone, extractAffyBatch, extractSnpFeatureSet, fromFiles, gcrmaSummary, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getData, getFullName, getIdentifier, getIntensities, getName, getSampleNames, getSiblings, getUnitIntensities, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, rmaSummary, setCdf, transformAffine

\bold{Methods inherited from AffymetrixFileSet}:\cr
append, appendFiles, as, as.AffymetrixFileSet, as.character, as.list, clearCache, clone, extract, fromFiles, getDescription, getFile, getFiles, getFileSize, getFullName, getIdentifier, getName, getNames, getPath, getPathnames, getTags, lapply, length, nbrOfFiles, seq, setName, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
