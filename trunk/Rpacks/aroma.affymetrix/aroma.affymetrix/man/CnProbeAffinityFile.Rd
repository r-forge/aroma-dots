%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnProbeAffinityFile}
\docType{class}
\alias{CnProbeAffinityFile}


 \title{The CnProbeAffinityFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnProbeAffinityFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeAffinityFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpProbeAffinityFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{CnProbeAffinityFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CnProbeAffinityFile}\cr
extends \link[aroma.affymetrix]{SnpProbeAffinityFile}\cr



  This class represents estimates of probe affinities in SNP probe-level models.
 }

 \usage{CnProbeAffinityFile(..., combineAlleles=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{SnpProbeAffinityFile}}.}
   \item{combineAlleles}{If \code{\link[base:logical]{FALSE}}, allele A and allele B are treated
      seperately, otherwise together.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{setCombineAlleles} \tab  -\cr
}


\bold{Methods inherited from SnpProbeAffinityFile}:\cr
getCellIndices, setMergeStrands

\bold{Methods inherited from ProbeAffinityFile}:\cr
getCellIndices, readUnits, updateUnits, writeSpatial

\bold{Methods inherited from ParameterCelFile}:\cr
decode, encode, readUnits, setDecodeFunction, setEncodeFunction, updateUnits

\bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, annotateMvsA, as.character, bgAdjustGcrma, bgAdjustOptical, bgAdjustRma, calibrateAllelicCrosstalk, clearData, clone, createFrom, fromFile, getAm, getCdf, getChipType, getData, getHeader, getHeaderV3, getIdentifier, getRectangle, getSampleName, getTimestamp, highlight, image270, importFeatureSet, nbrOfCells, normalizeQuantile, plotDensity, plotMvsA, plotMvsPosition, plotMvsX, range, readUnits, setCdf, smoothScatterMvsA, transformAffine, updateUnits, writeSpatial

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, clone, copyFile, fromFile, getFilename, getFileSize, getFileType, getFullName, getLabel, getName, getPath, getPathname, getTags, setLabel, stextLabel, stextLabels, stextSize

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
