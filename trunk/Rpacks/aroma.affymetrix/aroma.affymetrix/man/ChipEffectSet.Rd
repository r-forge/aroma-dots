%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ChipEffectSet}
\docType{class}
\alias{ChipEffectSet}


 \title{The ChipEffectSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ChipEffectSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{ChipEffectSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ChipEffectSet}\cr
extends \link[aroma.affymetrix]{AffymetrixCelSet}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{ChipEffectSet(..., model=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{readUnits} \tab  -\cr
}


\bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, applyToUnitIntensities, as, as.AffymetrixCelSet, as.character, fromFiles, getAverageFile, getCdf, getIntensities, getName, getUnitIntensities, nbrOfArrays, readUnits, setCdf

\bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, as.character, as.list, clearCache, extract, fromFiles, getFileSize, getName, getNames, getPath, lapply, length, nbrOfFiles, seq

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtained through the
   \code{getChipEffects()} method for the \code{\link{ProbeLevelModel}} class.
 }
\keyword{classes}
