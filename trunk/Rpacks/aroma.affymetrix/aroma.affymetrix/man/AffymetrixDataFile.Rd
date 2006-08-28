%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixDataFile}
\docType{class}
\alias{AffymetrixDataFile}


 \title{The abstract AffymetrixDataFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixDataFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{AffymetrixDataFile}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixApdFile}\cr

public abstract static class \bold{AffymetrixDataFile}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr



  An AffymetrixDataFile object represents a single Affymetrix data file,
  e.g. an Affymetrix CEL file or an Affymetrix Probe Data (APD) file.
  Note that this class is abstract and can not be instanciated, but
  instead you have to use one of the subclasses or the
  \code{\link[aroma.affymetrix:fromFile.AffymetrixDataFile]{*fromFile}()} method.
 }

 \usage{AffymetrixDataFile(...)}

 \arguments{
   \item{filename}{The filename of the data file.}
   \item{path}{An optional path to the data data file.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:000-LB-.AffymetrixDataFile]{[}} \tab Gets the probe intensities for a set of probes.\cr
 \tab \code{\link[aroma.affymetrix:as.character.AffymetrixDataFile]{as.character}} \tab Returns a short string describing the Affymetrix data file.\cr
 \tab \code{\link[aroma.affymetrix:fromFile.AffymetrixDataFile]{fromFile}} \tab Defines an AffymetrixDataFile object from a file.\cr
 \tab \code{\link[aroma.affymetrix:getProbeIntensities.AffymetrixDataFile]{getProbeIntensities}} \tab Gets the probe intensities for a set of probes.\cr
 \tab \code{\link[aroma.affymetrix:getProbesetIntensities.AffymetrixDataFile]{getProbesetIntensities}} \tab Gets probe signals for a subset of probesets in one array.\cr
 \tab \code{\link[aroma.affymetrix:normalizeQuantile.AffymetrixDataFile]{normalizeQuantile}} \tab Normalizes the sample to a target empirical distribution.\cr
 \tab \code{\link[aroma.affymetrix:plotDensity.AffymetrixDataFile]{plotDensity}} \tab Plots the densities of all samples.\cr
 \tab \code{transformAffine} \tab  -\cr
 \tab \code{writeApd} \tab  -\cr
}


\bold{Methods inherited from AffymetrixFile}:\cr
as.character, copyFile, fromFile, getFilename, getFileSize, getFileType, getLabel, getName, getPath, getPathname, setLabel, stextLabel, stextLabels

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
