%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffineProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffineProbeAffinityFile}
\docType{class}
\alias{AffineProbeAffinityFile}


 \title{The AffineProbeAffinityFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffineProbeAffinityFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{AffineProbeAffinityFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffineProbeAffinityFile}\cr
extends \emph{\link[aroma.affymetrix]{ProbeAffinityFile}}\cr



  This class represents estimates of probe affinities in the
  Bengtsson \& Hössjer (2006) model.
 }

 \usage{AffineProbeAffinityFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from ProbeAffinityFile}:\cr
findUnitsTodo, getFirstCellIndices, readUnits, updateUnits, writeSpatial

\bold{Methods inherited from ParameterCelFile}:\cr
decode, decodeUnit, decodeUnitGroup, encode, encodeUnit, encodeUnitGroup, readUnits, updateUnits

\bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, annotateMvsA, calcMvsA, clearData, clone, createFrom, fromFile, getCdf, getChipType, getData, getHeader, getRectangle, image270, nbrOfCells, normalizeQuantiles, plotMvsA, plotMvsPosition, plotMvsX, readUnits, setCdf, smoothScatterMvsA, stextChipType, transformAffine, updateUnits, writeSpatial

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, clone, copyFile, fromFile, getFilename, getFileSize, getFileType, getLabel, getName, getPath, getPathname, setLabel, stextLabel, stextLabels

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtain through the
   \code{getProbeAffinities()} method for the
   \code{\link{AffymetrixAffineModel}} class.
 }
\keyword{classes}
