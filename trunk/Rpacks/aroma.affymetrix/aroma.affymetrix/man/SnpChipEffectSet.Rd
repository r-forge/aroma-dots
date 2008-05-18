%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectSet.R
% 
% on Sat May 17 11:56:33 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SnpChipEffectSet}
\docType{class}
\alias{SnpChipEffectSet}


 \title{The SnpChipEffectSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SnpChipEffectSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{SnpChipEffectSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CnChipEffectSet}\cr

public static class \bold{SnpChipEffectSet}\cr
extends \link[aroma.affymetrix]{ChipEffectSet}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{SnpChipEffectSet(..., mergeStrands="byFirstFile")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ChipEffectSet}}.}
   \item{mergeStrands}{Specifies if the strands are merged or not for these
      estimates.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractTheta} \tab  -\cr
 \tab \code{extractTotalAndFreqB} \tab  -\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getBaseline} \tab  -\cr
 \tab \code{getMergeStrands} \tab  -\cr
 \tab \code{setMergeStrands} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectSet}:\cr
boxplotStats, calculateBaseline, calculateFieldBoxplotStats, calculateNuseBoxplotStats, calculateRleBoxplotStats, extractMatrix, extractTheta, findByName, findUnitsTodo, fromDataSet, fromFiles, getAM, getAsFullCelSet, getAverageFile, getBaseline, getCellIndices, getParameters, getXAM, plotBoxplot, readUnits, updateUnits

 \bold{Methods inherited from ParameterCelSet}:\cr
extractDataFrame, extractMatrix

 \bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, append, as, as.AffymetrixCelSet, averageQuantile, byName, extractAffyBatch, extractMatrix, findByName, fromFiles, fromName, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getChipType, getData, getFullName, getIntensities, getTimestamps, getUnitGroupCellMap, getUnitIntensities, isDuplicated, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, setCdf

 \bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, fromFiles, getPlatform

 \bold{Methods inherited from AromaMicroarrayDataSet}:\cr
getChipType, getPlatform, setAttributesBy, setAttributesBySampleAnnotationFile, setAttributesBySampleAnnotationSet, validate

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, equals, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getName, getNames, getPath, getPathnames, getTags, hasFile, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setTags, validate

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
