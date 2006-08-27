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
\cr

public static class \bold{ChipEffectFile}\cr
extends \emph{\link[aroma.affymetrix]{ParameterCelFile}}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{ChipEffectFile(..., model=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{createFrom} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


\bold{Methods inherited from ParameterCelFile}:\cr
createFrom, decode, decodeUnit, decodeUnitGroup, encode, encodeUnit, encodeUnitGroup, readUnits, updateUnits

\bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, fromFile, getCdf, getChipType, getFields, getHeader, getRectangle, image270, normalizeQuantiles, readUnits, setCdf, transformAffine, updateUnits, writeSpatial

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, fromFile, getFilename, getFileSize, getFileType, getName, getPath, getPathname

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtained through the
   \code{getChipEffects()} method for the \code{\link{ProbeLevelModel}} class.
 }
\keyword{classes}
