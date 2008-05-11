%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ExonChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ExonChipEffectFile}
\docType{class}
\alias{ExonChipEffectFile}


 \title{The ExonChipEffectFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ExonChipEffectFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{ExonChipEffectFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ExonChipEffectFile}\cr
extends \link[aroma.affymetrix]{ChipEffectFile}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{ExonChipEffectFile(..., mergeGroups=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ChipEffectFile}}.}
   \item{mergeGroups}{Specifies if the groups are merged or not for these
      estimates.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectFile}:\cr
extractMatrix, extractTheta, findUnitsTodo, getAM, getAsFullCelFile, getCellIndices, getCellMapForMainCdf, getExpandedCellMap, getParameters, getUnitGroupCellChromosomePositionMap, getXAM, mergeGroups, readUnits, writeAsFullCelFile

 \bold{Methods inherited from ParameterCelFile}:\cr
extractDataFrame, extractMatrix, getParameters, readUnits

 \bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, allocateFromCdf, clone, createFrom, fromFile, getAm, getAttributeXY, getCdf, getFileFormat, getImage, getPloidy, hasAttributeXY, highlight, image270, nbrOfCells, plotDensity, plotImage, plotMvsA, plotMvsX, setAttributesByTags, setAttributeXY, setCdf, smoothScatterMvsA, writeImage

 \bold{Methods inherited from AffymetrixFile}:\cr
getPlatform

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getChipType, getPlatform

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of a \code{\link{ExonChipEffectSet}}.
 }
\keyword{classes}
