%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawData.R
% 
% on Tue Jan 15 18:35:57 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{RawData}
\docType{class}
\alias{RawData}


\title{The RawData class}

\description{
 Package:  aroma \cr
\bold{Class RawData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{RawData}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{RawData}\cr
extends \link[aroma]{MicroarrayData}\cr



 Creates a new \code{RawData} object.
}

\usage{RawData(R=NULL, G=NULL, Rb=NULL, Gb=NULL, layout=NULL, extras=list())}

\arguments{
  \item{R,G}{A NxM \code{\link[base]{matrix}} containing (non-logged) foreground signals of
   the red (green) channel, where N is the number of spots on each slide
   and M is the number of slides in this data set.}
  \item{Rb,Gb}{A NxM \code{\link[base]{matrix}} containing (non-logged) background signals of
   the red (green) channel, where N is the number of spots on each slide
   and M is the number of slides in this data set.}
  \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
   slides in this data set.}
  \item{extras}{Private argument. Do not use.}
}

\section{Fields and Methods}{
 \tabular{rll}{
   \tab \code{R} \tab The foreground for channel R (non-logged). \cr
   \tab \code{G} \tab The foreground for channel G (non-logged). \cr
   \tab \code{Rb} \tab The background for channel Rb (non-logged). \cr
   \tab \code{Gb} \tab The background for channel Gb (non-logged). \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.RawData} \tab  -\cr
 \tab \code{as.RGData} \tab  -\cr
 \tab \code{boxplot} \tab  -\cr
 \tab \code{calibrateMultiscan} \tab  -\cr
 \tab \code{\link[aroma:getBackground.RawData]{getBackground}} \tab Gets the background signal.\cr
 \tab \code{\link[aroma:getColors.RawData]{getColors}} \tab Generates colors for each of the specified spots.\cr
 \tab \code{\link[aroma:getForeground.RawData]{getForeground}} \tab Gets the foreground signal.\cr
 \tab \code{\link[aroma:getSignal.RawData]{getSignal}} \tab Gets the red and the green signal.\cr
 \tab \code{getSNR1} \tab  -\cr
 \tab \code{getSNR2} \tab  -\cr
 \tab \code{getWithinChannelPairs} \tab  -\cr
 \tab \code{normalizeAffine} \tab  -\cr
 \tab \code{normalizeGenewise} \tab  -\cr
 \tab \code{normalizeQuantile} \tab  -\cr
 \tab \code{normalizeWithinSlide} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{plotSpatial} \tab  -\cr
 \tab \code{plotXY} \tab  -\cr
 \tab \code{range} \tab  -\cr
 \tab \code{\link[aroma:shift.RawData]{shift}} \tab Shift the log-ratios, log-intensities or the raw signal.\cr
 \tab \code{\link[aroma:swapDyes.RawData]{swapDyes}} \tab Dye swap one or many slides.\cr
}


 \bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\note{
 There are several functions that returns an object of this class, and
 it is only in very special cases that you actually have to create one
 yourself.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{For an example see help(MAData).}}
\keyword{classes}
