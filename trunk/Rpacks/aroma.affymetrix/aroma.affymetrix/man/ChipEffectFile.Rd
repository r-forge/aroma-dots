%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ChipEffectFile}
\docType{class}
\alias{ChipEffectFile}


 \title{The ChipEffectFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ChipEffectFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{ChipEffectFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CnChipEffectFile}, \link[aroma.affymetrix]{ExonChipEffectFile}, \link[aroma.affymetrix]{SnpChipEffectFile}\cr

public static class \bold{ChipEffectFile}\cr
extends \link[aroma.affymetrix]{ParameterCelFile}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{ChipEffectFile(..., probeModel=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:extractDataFrame.ChipEffectFile]{extractDataFrame}} \tab Extract data as a data frame.\cr
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getAM.ChipEffectFile]{getAM}} \tab Gets the log-intensities and log-ratios of chip effects for two arrays.\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getXAM.ChipEffectFile]{getXAM}} \tab Gets the physical position, log-intensities and log-ratios of chip effects for two arrays.\cr
 \tab \code{readUnits} \tab  -\cr
}


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

 \seealso{
   An object of this class is typically obtained through the
   \code{getChipEffects()} method for the \code{\link{ProbeLevelModel}} class.
   An object of this class is typically part of a \code{\link{ChipEffectSet}}.
 }
\keyword{classes}
