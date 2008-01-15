%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyze20Data.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ScanAlyze20Data}
\docType{class}
\alias{ScanAlyze20Data}


\title{The ScanAlyze20Data (ScanAlyze v2.0) class}

\usage{ScanAlyze20Data(layout=NULL)}

\description{
 Package:  aroma \cr
\bold{Class ScanAlyze20Data}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{ScanAlyzeData}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{ScanAlyze20Data}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ScanAlyze20Data}\cr
extends \link[aroma]{ScanAlyzeData}\cr



 Creates an empty ScanAlyze20Data object.
}

\arguments{
  \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
   slides in this data set.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getArea} \tab  -\cr
 \tab \code{getBackground} \tab  -\cr
 \tab \code{getBackgroundSD} \tab  -\cr
 \tab \code{getBgArea} \tab  -\cr
 \tab \code{getDiameter} \tab  -\cr
 \tab \code{getForeground} \tab  -\cr
 \tab \code{getForegroundSD} \tab  -\cr
 \tab \code{getGrid} \tab  -\cr
 \tab \code{getSpotColumn} \tab  -\cr
 \tab \code{\link[aroma:getSpotPosition.ScanAlyze20Data]{getSpotPosition}} \tab Gets physical positions of the spots.\cr
 \tab \code{getSpotRow} \tab  -\cr
 \tab \code{\link[aroma:read.ScanAlyze20Data]{read}} \tab Reads one or several ScanAlyze v2.0 files.\cr
 \tab \code{\link[aroma:write.ScanAlyze20Data]{write}} \tab Writes a ScanAlyze v2.0 Results Data file (NOT IMPLEMENTED).\cr
}


\bold{Methods inherited from ScanAlyzeData}:\cr
append, as.RawData, getArea, getBackground, getBgArea, getCircularity, getDiameter, getForeground, getGrid, getRawData, getSpotColumn, getSpotPosition, getSpotRow, normalizeGenewise, plotSpatial, plotSpatial3d, read, write

\bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] \url{http://array.mbb.yale.edu/analysis/manual/inputfiles/}
}
\keyword{classes}
