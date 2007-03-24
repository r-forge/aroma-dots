%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnChipEffectFile}
\docType{class}
\alias{CnChipEffectFile}


 \title{The CnChipEffectFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnChipEffectFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpChipEffectFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{CnChipEffectFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CnChipEffectFile}\cr
extends \link[aroma.affymetrix]{SnpChipEffectFile}\cr



  This class represents estimates of chip effects in a copy-number probe-level
  models.
 }

 \usage{CnChipEffectFile(..., combineAlleles=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{SnpChipEffectFile}}.}
   \item{combineAlleles}{A \code{\link[base]{logical}} indicating if the signals from allele A
     and allele B are combined or not.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{mergeStrands} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


 \bold{Methods inherited from SnpChipEffectFile}:\cr
getCellIndices, getParameters, mergeStrands, readUnits

 \bold{Methods inherited from ChipEffectFile}:\cr
extractDataFrame, findUnitsTodo, getAM, getCellIndices, getParameters, getXAM, mergeGroups, readUnits

 \bold{Methods inherited from ParameterCelFile}:\cr
readUnits

 \bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, clone, createFrom, fromFile, getAm, getAttributeXY, getCdf, getImage, getPloidy, hasAttributeXY, highlight, image270, nbrOfCells, plotDensity, plotImage, plotMvsA, plotMvsPosition, plotMvsX, setAttributesByTags, setAttributeXY, setCdf, smoothScatterMvsA, writeImage

 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of a \code{\link{CnChipEffectSet}}.
 }
\keyword{classes}
