%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{RGData}
\docType{class}
\alias{RGData}


\title{The RGData class}

\description{
 Package:  aroma \cr
\bold{Class RGData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{RGData}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{RGData}\cr
extends \link[aroma]{MicroarrayData}\cr



 Creates a new \code{RGData} object.
 The philosophy behind this data structure is to think about the data in the form of the signals in one channel (R) versus the signals in the other channel (G).
 This is in contrast to the idea of the \code{\link{MAData}} structure, which thinks about the data as the log ratios (M) and log intensites (A) for the spot signals.
}

\usage{RGData(R=NULL, G=NULL, layout=NULL, extras=list())}

\arguments{
  \item{R,G}{A NxM \code{\link[base]{matrix}} containing (non-logged) signals of the red
   (green) channel, where N is the number of spots on each slide and M
   is the number of slides in this data set.}
  \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
   slides in this data set.}
  \item{extras}{Private argument. Do not use.}
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
   \tab \code{R} \tab The signal for channel R (non-logged). \cr
   \tab \code{G} \tab The signal for channel G (non-logged). \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma:as.MAData.RGData]{as.MAData}} \tab Transform from the red and green intensities into log ratios between them and the log product of them.\cr
 \tab \code{as.RawData} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.RGData} \tab  -\cr
 \tab \code{boxplot} \tab  -\cr
 \tab \code{\link[aroma:calibrateMultiscan.RGData]{calibrateMultiscan}} \tab Calibrates multiple re-scanned images based on an affine model.\cr
 \tab \code{\link[aroma:fitMultiscanAffine.RGData]{fitMultiscanAffine}} \tab Fits an affine model to signals from a multi-scanned slide.\cr
 \tab \code{getA} \tab  -\cr
 \tab \code{getCalibratedMultiscan} \tab  -\cr
 \tab \code{getChannelNames} \tab  -\cr
 \tab \code{\link[aroma:getColors.RGData]{getColors}} \tab Generates red to green colors for each of the specified spots.\cr
 \tab \code{\link[aroma:getLogIntensities.RGData]{getLogIntensities}} \tab Calculates the log-intensitites (A values).\cr
 \tab \code{\link[aroma:getLogRatios.RGData]{getLogRatios}} \tab Calculates the log-ratios (M values).\cr
 \tab \code{getM} \tab  -\cr
 \tab \code{getWithinChannelPairs} \tab  -\cr
 \tab \code{log} \tab  -\cr
 \tab \code{\link[aroma:mean.RGData]{mean}} \tab Genewise Average Mean for channel R and G.\cr
 \tab \code{\link[aroma:normalizeAffine.RGData]{normalizeAffine}} \tab Weighted affine normalization between channels and arrays.\cr
 \tab \code{\link[aroma:normalizeCurveFit.RGData]{normalizeCurveFit}} \tab Within-slide intensity-dependent normalization in (A,M).\cr
 \tab \code{normalizeGenewise} \tab  -\cr
 \tab \code{normalizeLoess} \tab  -\cr
 \tab \code{normalizeLowess} \tab  -\cr
 \tab \code{normalizeQuantile} \tab  -\cr
 \tab \code{normalizeRobustSpline} \tab  -\cr
 \tab \code{normalizeSpline} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{plotSpatial} \tab  -\cr
 \tab \code{plotXY} \tab  -\cr
 \tab \code{power} \tab  -\cr
 \tab \code{range} \tab  -\cr
 \tab \code{\link[aroma:shift.RGData]{shift}} \tab Shift the log-ratios, log-intensities or the raw signal.\cr
 \tab \code{\link[aroma:swapDyes.RGData]{swapDyes}} \tab Swap dyes of one or many slides.\cr
 \tab \code{\link[aroma:var.RGData]{var}} \tab Genewise Variance for channel R and G.\cr
}


\bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}

\details{
  The mapping between M and A, and R and G is a one-to-one function.
  Given the signal R and G for the R and the G channels you get the
  M and the A values by:
  \deqn{
    M = \log_2\frac{R}{G},\quad
    A = \log_2\sqrt{R{\cdot}G} = \frac{1}{2}\log_2 R{\cdot}G,
  }{
    M = log2(R/G), A = log2(sqrt(R*G)) = 1/2*log2(R*G),
  }
  which in [R] can be done by \code{ma <- as.MAData(rg)}. The reverse function, i.e. going back to the R and the G is:
  \deqn{
    R = \sqrt{2^{2A+M}},\quad G = \sqrt{2^{2A-M}}
  }{
    R = sqrt(2^(2A+M)), G = sqrt(2^(2A-M))
  }
  which in [R] can be done by \code{rg <- as.RGData(rg)}.

  Note that if the signal in one or both channels is non-positive,
  the log-transform will make these values undefined, that is, set
  them to \code{\link[base]{NA}}. When going back to (G,R) from (A,M) these values
  will remain \code{\link[base]{NA}}.
}

\note{
 There are several functions that returns an object of this class, and
 it is only in very special cases that you actually have to create one
 yourself.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  # Create a raw data object from the preexisting example data in
  # the sma package.
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)

  # Get the signal (here by default non-background corrected)
  ma <- getSignal(raw)

  # Transform (M,A) into (R,G)
  rg <- as.RGData(ma)
}
\keyword{classes}
