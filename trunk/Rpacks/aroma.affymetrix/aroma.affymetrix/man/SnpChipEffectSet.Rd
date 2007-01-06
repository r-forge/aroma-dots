%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectSet.R
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
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectSet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{SnpChipEffectSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CnChipEffectSet}\cr

public static class \bold{SnpChipEffectSet}\cr
extends \link[aroma.affymetrix]{ChipEffectSet}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{SnpChipEffectSet(..., mergeStrands=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ChipEffectSet}}.}
   \item{mergeStrands}{Specifies if the strands are merged or not for these
      estimates.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:extractSnpQSet.SnpChipEffectSet]{extractSnpQSet}} \tab Extracts chip effects as a SnpQSet.\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getMergeStrands} \tab  -\cr
 \tab \code{setMergeStrands} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectSet}:\cr
findUnitsTodo, fromDataSet, fromFiles, getAverageFile, getCellIndices, readUnits, updateUnits

 \bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, append, as, as.AffymetrixCelSet, averageQuantile, cna, extractAffyBatch, extractSnpFeatureSet, fromFiles, gcrmaSummary, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getData, getFullName, getIntensities, getName, getTimestamps, getUnitIntensities, isDuplicated, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, rmaSummary, setCdf

 \bold{Methods inherited from AffymetrixFileSet}:\cr
append, appendFiles, as, as.AffymetrixFileSet, as.list, clearCache, extract, fromFiles, getFile, getFileSize, getFullName, getName, getNames, getPath, getPathnames, getTags, lapply, nbrOfFiles, seq, setName, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
