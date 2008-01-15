%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{MAData}
\docType{class}
\alias{MAData}


\title{The MAData class}

\description{
 Package:  aroma \cr
\bold{Class MAData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{MAData}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{BMAData}, \link[aroma]{TMAData}\cr

public static class \bold{MAData}\cr
extends \link[aroma]{MicroarrayData}\cr



 Creates a new \code{MAData} object.
 The philosophy behind this data structure is to think about the data in the form of log ratios (M) and log intensites (A) for the spot signals.
 This is in contrast to the idea of the \link{RGData} structure, which thinks about the data as the signal in one channel (R) versus the signal in the other channel (G).
}

\usage{MAData(M=NULL, A=NULL, layout=NULL, extras=list())}

\arguments{
  \item{M,A}{A NxM \code{\link[base]{matrix}} containing (base-2) log-ratios and
   log-intensities, respectively.}
  \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
   slides in this data set.}
  \item{extras}{Private argument. Do not use.}
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
   \tab \code{M} \tab The log of the ratios between R and G, i.e. \eqn{\log_2{\frac{R}{G}}}{log2(R/G)} where R and G is the signal for the R and the G channel. \cr
   \tab \code{A} \tab The log of the intensities of R and G, i.e. \eqn{\frac{1}{2}\log_2{R{\cdot}G}}{1/2*log2(R*G)} where R and G is the signal for the R and the G channel. \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.MAData} \tab  -\cr
 \tab \code{as.RawData} \tab  -\cr
 \tab \code{\link[aroma:as.RGData.MAData]{as.RGData}} \tab Transform MA format into RG format.\cr
 \tab \code{boxplot} \tab  -\cr
 \tab \code{drawCurveFit} \tab  -\cr
 \tab \code{\link[aroma:estimateMACurve.MAData]{estimateMACurve}} \tab Estimates a smooth intensity-dependent curve in (A,M).\cr
 \tab \code{\link[aroma:getAdjustedSpotVariability.MAData]{getAdjustedSpotVariability}} \tab Gets the spotwise intensity-adjusted variability of replicate slides.\cr
 \tab \code{\link[aroma:getColors.MAData]{getColors}} \tab Generates colors for each of the specified spots.\cr
 \tab \code{getDenseSpots} \tab  -\cr
 \tab \code{getG} \tab  -\cr
 \tab \code{\link[aroma:getGeneVariability.MAData]{getGeneVariability}} \tab Gets the genewise variability.\cr
 \tab \code{getHistogram} \tab  -\cr
 \tab \code{\link[aroma:getMOR.MAData]{getMOR}} \tab Gets the Measure of Reproducibility.\cr
 \tab \code{\link[aroma:getMOR2003a.MAData]{getMOR2003a}} \tab Gets the Intensity Adjusted Measure of Reproducibility.\cr
 \tab \code{getR} \tab  -\cr
 \tab \code{getRange} \tab  -\cr
 \tab \code{getSpotVariability} \tab  -\cr
 \tab \code{hist} \tab  -\cr
 \tab \code{loessQuantile} \tab  -\cr
 \tab \code{\link[aroma:mean.MAData]{mean}} \tab Average Mean for microarray data.\cr
 \tab \code{\link[aroma:normalizeAcrossSlides.MAData]{normalizeAcrossSlides}} \tab Normalizes across slides.\cr
 \tab \code{\link[aroma:normalizeAffine.MAData]{normalizeAffine}} \tab Affine normalization based on non-logged data.\cr
 \tab \code{\link[aroma:normalizeCurveFit.MAData]{normalizeCurveFit}} \tab Within-slide normalization that adjust log-ratios by estimating a smooth intensity-dependent curve in (A,M).\cr
 \tab \code{normalizeGenewise} \tab  -\cr
 \tab \code{\link[aroma:normalizeLogRatioShift.MAData]{normalizeLogRatioShift}} \tab Within-slide normalization that adjust the log-ratios shift.\cr
 \tab \code{normalizePlatewise} \tab  -\cr
 \tab \code{normalizePrintorder} \tab  -\cr
 \tab \code{normalizeQuantile} \tab  -\cr
 \tab \code{normalizeSpatial} \tab  -\cr
 \tab \code{\link[aroma:normalizeWithinSlide.MAData]{normalizeWithinSlide}} \tab Within-slide normalization.\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{plot3d} \tab  -\cr
 \tab \code{plotDensity} \tab  -\cr
 \tab \code{plotDiporder} \tab  -\cr
 \tab \code{\link[aroma:plotMvsM.MAData]{plotMvsM}} \tab Plots the log-ratios for one slide versus another.\cr
 \tab \code{plotPrintorder} \tab  -\cr
 \tab \code{plotSpatial} \tab  -\cr
 \tab \code{plotSpatial3d} \tab  -\cr
 \tab \code{plotXY} \tab  -\cr
 \tab \code{qqnorm} \tab  -\cr
 \tab \code{range} \tab  -\cr
 \tab \code{\link[aroma:shift.MAData]{shift}} \tab Shift the log-ratios, log-intensities or the raw signal.\cr
 \tab \code{shiftEqualRG} \tab  -\cr
 \tab \code{\link[aroma:swapDyes.MAData]{swapDyes}} \tab Dye swap one or many slides.\cr
 \tab \code{\link[aroma:topSpots.MAData]{topSpots}} \tab Gets the top spots.\cr
 \tab \code{\link[aroma:var.MAData]{var}} \tab Variance for microarray data.\cr
}


\bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}

\details{
  The mapping between M and A, and R and G is a one-to-one function, i.e.
  you can go back and forth without loosing any information.
  Given the signal R and G for the R and the G channels you get the
  M and the A values by:
  \deqn{
    M = \log_2\frac{R}{G},\quad
    A = \log_2\sqrt{R{\cdot}G} = \frac{1}{2}\log_2 R{\cdot}G,
  }{
    M = log2(R/G), A = log2(sqrt(R*G)) = 1/2*log2(R*G),
  }
  and going back to the R and the G by:
  \deqn{
    R = \sqrt{2^{2A+M}},\quad G = \sqrt{2^{2A-M}}
  }{
    R = sqrt(2^(2A+M)), G = sqrt(2^(2A-M))
  }

  Comments: M is a memnonic for Minus and A is for Add since
  \eqn{M = \log_2{R}-\log_2{G}} and \eqn{A = 1/2*(\log_2{R}+\log_2{G})}.
}

\note{
 There are several functions that returns an object of this class, and
 it is only in very special cases that you actually have to create one
 yourself.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:mouse.data")

  # Create a raw data object from the preexisting example data in
  # the sma package.
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)

  # Get the signal (here by default non-background corrected)
  ma <- getSignal(raw, bgSubtract=TRUE)

  # Transform (M,A) into (R,G)
  rg <- as.RGData(ma)

  # Transform back from (R,G) to (M,A)
  ma2 <- as.MAData(rg);

  # Check that the tranformation a one-to-one function
  print(equals(ma, ma2))   # TRUE

  layout(matrix(1:4, ncol=2, byrow=TRUE))
  # Plot the R vs G with a fitted (lowess) line for slide 2.
  plot(rg, slide=2); lowessCurve(rg)

  # And the similar for M vs A.
  plot(ma, slide=2); lowessCurve(ma)

  # Plot a spatial representation of the M's.
  plotSpatial(ma, slide=2)

  # Make a boxplot of the print-tip groups.
  boxplot(ma, groupBy="printtip", slide=2)
}
\keyword{classes}
