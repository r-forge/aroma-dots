%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterCelSet.R
% 
% on Sat May 17 11:55:25 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ParameterCelSet}
\docType{class}
\alias{ParameterCelSet}


 \title{The ParameterCelSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ParameterCelSet}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{ParameterCelSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ChipEffectSet}, \link[aroma.affymetrix]{CnChipEffectSet}, \link[aroma.affymetrix]{ExonChipEffectSet}, \link[aroma.affymetrix]{FirmaSet}, \link[aroma.affymetrix]{SnpChipEffectSet}\cr

public static class \bold{ParameterCelSet}\cr
extends \link[aroma.affymetrix]{AffymetrixCelSet}\cr



  A ParameterCelSet object represents a set of \code{\link{ParameterCelFile}}:s.
 }

 \usage{ParameterCelSet(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixCelSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractDataFrame} \tab  -\cr
 \tab \code{extractMatrix} \tab  -\cr
}


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
\keyword{IO}
