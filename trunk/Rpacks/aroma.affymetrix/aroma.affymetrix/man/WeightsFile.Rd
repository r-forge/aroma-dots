%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  WeightsFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{WeightsFile}
\docType{class}
\alias{WeightsFile}


 \title{The WeightsFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class WeightsFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{WeightsFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{WeightsFile}\cr
extends \link[aroma.affymetrix]{ParameterCelFile}\cr



  This class represents weights calculated from residuals of
  probe-level models.
 }

 \usage{WeightsFile(..., probeModel=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getImage} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{writeImage} \tab  -\cr
}


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
   An object of this class is typically obtained through the
   \code{getWeightsSet()} method for the \code{\link{ProbeLevelModel}} class.
   An object of this class is typically part of a \code{\link{WeightsSet}}.
 }
\keyword{classes}
