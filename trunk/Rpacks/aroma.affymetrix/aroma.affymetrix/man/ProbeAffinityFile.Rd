%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ProbeAffinityFile}
\docType{class}
\alias{ProbeAffinityFile}


 \title{The ProbeAffinityFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ProbeAffinityFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{ProbeAffinityFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{SnpProbeAffinityFile}\cr

public static class \bold{ProbeAffinityFile}\cr
extends \link[aroma.affymetrix]{ParameterCelFile}\cr



  This class represents estimates of probe affinities in probe-level models.
 }

 \usage{ProbeAffinityFile(..., model=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
   \item{model}{The specific type of model, e.g. \code{"pm"}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{writeSpatial} \tab  -\cr
}


\bold{Methods inherited from ParameterCelFile}:\cr
decode, encode, readUnits, setDecodeFunction, setEncodeFunction, updateUnits

\bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, annotateMvsA, bgAdjustGcrma, bgAdjustOptical, bgAdjustRma, calibrateAllelicCrosstalk, clearData, clone, createFrom, fromFile, getAm, getCdf, getChipType, getData, getHeader, getIdentifier, getRectangle, getSampleName, highlight, image270, nbrOfCells, normalizeQuantile, plotDensity, plotMvsA, plotMvsPosition, plotMvsX, range, readUnits, setCdf, smoothScatterMvsA, transformAffine, updateUnits, writeSpatial

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, clone, copyFile, fromFile, getFilename, getFileSize, getFileType, getFullName, getLabel, getName, getPath, getPathname, getTags, setLabel, stextLabel, stextLabels, stextSize

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtained through the
   \code{getProbeAffinities()} method for the \code{\link{ProbeLevelModel}} class.
 }
\keyword{classes}
