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
getCellIndices, getParameters, readUnits

 \bold{Methods inherited from ParameterCelFile}:\cr
readUnits

 \bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, clone, fromFile, getAm, getAttributeXY, getCdf, getImage, getPloidy, hasAttributeXY, highlight, image270, nbrOfCells, plotDensity, plotImage, plotMvsA, plotMvsPosition, plotMvsX, setAttributesByTags, setAttributeXY, setCdf, smoothScatterMvsA, writeImage

 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, fromFile, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, setAttribute, setAttributes, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
