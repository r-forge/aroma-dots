%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SnpChipEffectFile}
\docType{class}
\alias{SnpChipEffectFile}


 \title{The SnpChipEffectFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SnpChipEffectFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{SnpChipEffectFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CnChipEffectFile}\cr

public static class \bold{SnpChipEffectFile}\cr
extends \link[aroma.affymetrix]{ChipEffectFile}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{SnpChipEffectFile(..., mergeStrands=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ChipEffectFile}}.}
   \item{mergeStrands}{Specifies if the strands are merged or not for these
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
findUnitsTodo, getAM, getCellIndices, getParameters, getXAM, mergeGroups, readUnits

 \bold{Methods inherited from ParameterCelFile}:\cr
readUnits

 \bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, clone, fromFile, getAm, getCdf, getImage, highlight, image270, nbrOfCells, plotDensity, plotImage, plotMvsA, plotMvsPosition, plotMvsX, setCdf, smoothScatterMvsA, writeImage

 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, fromFile, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of a \code{\link{SnpChipEffectSet}}.
 }
\keyword{classes}
