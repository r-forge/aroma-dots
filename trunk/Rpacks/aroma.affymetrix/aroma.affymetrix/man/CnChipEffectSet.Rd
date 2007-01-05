%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnChipEffectSet.R
% 
% on Sat Jan 06 01:50:14 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnChipEffectSet}
\docType{class}
\alias{CnChipEffectSet}


 \title{The CnChipEffectSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnChipEffectSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectSet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpChipEffectSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{CnChipEffectSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CnChipEffectSet}\cr
extends \link[aroma.affymetrix]{SnpChipEffectSet}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{CnChipEffectSet(..., combineAlleles=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{SnpChipEffectSet}}.}
   \item{combineAlleles}{A \code{\link[base]{logical}} indicating if the signals from
      allele A and allele B are combined or not.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getCombineAlleles} \tab  -\cr
 \tab \code{setCombineAlleles} \tab  -\cr
 \tab \code{writeWig} \tab  -\cr
}


\bold{Methods inherited from SnpChipEffectSet}:\cr
extractSnpQSet, getAverageFile, getChipEffectFileClass, getMergeStrands, setMergeStrands

\bold{Methods inherited from ChipEffectSet}:\cr
findUnitsTodo, fromDataSet, fromFiles, getAverageFile, getCellIndices, getChipEffectFileClass, readUnits, updateUnits

\bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, append, applyToUnitIntensities, as.AffymetrixCelSet, as.character, averageQuantile, bgAdjustGcrma, bgAdjustOptical, bgAdjustRma, calculateGsbParameters, calibrateAllelicCrosstalk, clone, cna, createBlankSet, createFromFeatureSet, extractAffyBatch, extractSnpFeatureSet, fromFiles, gcrmaSummary, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getData, getFullName, getIdentifier, getIntensities, getName, getSiblings, getTimestamps, getUnitIntensities, isDuplicated, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, rmaSummary, setCdf, transformAffine

\bold{Methods inherited from AffymetrixFileSet}:\cr
append, appendFiles, as.AffymetrixFileSet, as.character, as.list, clearCache, clone, extract, fromFiles, getDescription, getFile, getFiles, getFileSize, getFullName, getIdentifier, getName, getNames, getPath, getPathnames, getTags, lapply, length, nbrOfFiles, reorder, seq, setName, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
\keyword{classes}
