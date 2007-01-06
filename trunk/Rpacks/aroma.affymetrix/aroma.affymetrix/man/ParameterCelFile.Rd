%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ParameterCelFile}
\docType{class}
\alias{ParameterCelFile}


 \title{The ParameterCelFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ParameterCelFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{ParameterCelFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ChipEffectFile}, \link[aroma.affymetrix]{CnChipEffectFile}, \link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{ProbeAffinityFile}, \link[aroma.affymetrix]{SnpChipEffectFile}, \link[aroma.affymetrix]{SnpProbeAffinityFile}\cr

public static class \bold{ParameterCelFile}\cr
extends \link[aroma.affymetrix]{AffymetrixCelFile}\cr



  A ParameterCelFile object represents parameter estimates.
 }

 \usage{ParameterCelFile(..., encodeFunction=NULL, decodeFunction=NULL)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
   \item{encodeFunction}{A \code{\link[base]{function}} taking a single \code{\link[base]{list}} structure
      as its argument.}
   \item{decodeFunction}{A \code{\link[base]{function}} taking a single \code{\link[base]{list}} structure
      as its argument.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{readUnits} \tab  -\cr
}


 \bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, clone, fromFile, getAm, getCdf, highlight, image270, nbrOfCells, plotDensity, plotMvsA, plotMvsPosition, plotMvsX, setCdf, smoothScatterMvsA

 \bold{Methods inherited from AffymetrixFile}:\cr
fromFile, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{File format}{
   The idea behind this class is store data fields which by nature have
   one value per probe (per field) in CEL files.  A perfect example is to
   store probe-affinity estimates and their standard deviations.  There
   is one probe affinity per probe so the structure of a CEL file (and
   its coupled CDF file) is well suited to read/write such information.
 
   Consider a unit group with L probes.  A CEL file stores
   \code{intensities} (L floats), \code{stdvs} (L floats), and
   \code{pixels} (L integers).  Thus, for each probe l=1,...,L, a
   (float, float, integer) tuple is stored.  We can use this for any
   information we want.  If we want a slightly different structure,
   we can choose to encode/decode our structure/information to fit the
   structure of the CEL file.  This abstract class provides transparent
   methods for encoding and decoding such information through methods
   \code{\link[aroma.affymetrix:encodeUnitGroup.ParameterCelFile]{*encodeUnitGroup}()} and \code{\link[aroma.affymetrix:decodeUnitGroup.ParameterCelFile]{*decodeUnitGroup}()}.
   By subclassing you can implement different types of data structures.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{classes}
\keyword{IO}
