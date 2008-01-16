%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantArrayData.R
% 
% on Tue Jan 15 18:35:00 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QuantArrayData}
\docType{class}
\alias{QuantArrayData}


\title{The QuantArrayData class}

\usage{QuantArrayData(layout=NULL)}

\description{
 Package:  aroma \cr
\bold{Class QuantArrayData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{QuantArrayData}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{QuantArrayData}\cr
extends \link[aroma]{MicroarrayData}\cr



 Creates an empty QuantArrayData object.
}

\arguments{
  \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
   slides in this data set.}
}

\section{Fields}{
  \item{ch1 Intensity, ch2 Intensity, ...}{The spot intensity.
    What the authors mean by intensity is not clear, but on p46 in [2]
    it says "Quantitation. Output allows the user to choose the way
    the spot intensities are calculated: a) the sum of the intensities
    for the various pixels of the spot and the background (Total
    Intensities), b) the average of the intensities (Mean Intensity),
    c) the most frequent occurrence of a pixel value (Mode Intensity),
    or d) the Median Intensity can be reported."

    If this is true, then where is the information about which method
    was used? Looking at the header of a typical QuantArray data file
    this information is not stored anywhere.
  }
  \item{ch1 Intensity Std Dev, ch2 Intensity Std Dev, ...}{
    The standard deviation of the spot pixel intensities.}
  \item{ch1 Spot Uniformity, ch2 Spot Uniformity, ...}{
    Value in [0,1].
    The 0.99 quantile minus the 0.01 quantile pixel signal divided
    by the pixel intensity range.}

  \item{ch1 Background, ch2 Background, ...}{The background intensity.}
  \item{ch1 Background Std Dev, ch2 Background Std Dev, ...}{
    The standard deviation of the background pixel intensities.}
  \item{ch1 Bkg. Uniformity, ch2 Bkg. Uniformity, ...}{Value in [0,1].}

  \item{ch1 Diameter, ch2 Diameter, ...}{The average spot diameter
    (in micrometers).}
  \item{ch1 Area, ch2 Area, ...}{The spot area as the number of pixels
    (in square micrometers).}
  \item{ch1 Footprint, ch2 Footprint, ...}{The distance (in micrometers)
    between the expected position of the spot and its actual position.}
  \item{ch1 Circularity, ch2 Circularity, ...}{
    Value in [0,1]. The circularity.}
  \item{ch1 Replicate Uniformity, ch2 Replicate Uniformity, ...}{
    Value in [0,1]. Uniformity of replicated spots.}

  \item{ch1 Confident, ch2 Confident, ...}{
    Value in [0,1].
    Different methods may have been used. See in
    \code{sections[["Protocol Info"]]} of the QuantArrayData object
    which was used.
    Let \code{(c[1], c[2], ..., c[Q])} be the individual measurement
    confidences. Then
    "Minimum") \code{cAll = min(c)}.
    "Weighted Average") \code{cAll = sum(w*c)/sum(w)} where
    \code{w} are weights.
    "Product") \code{cAll = (prod(c))^(1/Q)}.
  }

 \emph{Redundant fields}:\cr
  \item{ch1 Percentage, ch2 Percentage, ...}{\bold{REDUNDANT!}
    \code{=="ch1 Intensity"/("ch1 Intensity"+"ch2 Intensity")} and so on.
    For each spot, the percentage of the total brightness (the sum of all
    channels' intensities) for each channel.}
  \item{ch1 Ratio, ch2 Ratio, ...}{\bold{REDUNDANT!}
    \code{=="ch2 Intensity"/"ch1 Intensity"} or similar.
    For each spot, the relative intensity of a spot, compared to the
    corresponding spot in the control. The ratio value for the channel
    that is the control channel (commonly channel 1) are all equal to 1.
    The ratio value for the other channels are equal to their relative
    ratio of the intensity compared to the control.
    Background subtraction might have been applied first, but that is
    still reduandant.}
  \item{ch1 Signal Noise Ratio, ch2 Signal Noise Ratio, ...}{
    \bold{REDUNDANT!}
    \code{=="ch1 Intensity"/"ch1 Background Std Dev"} and so on.
    The ratio of the spot intensity to the standard deviation of the
    local background of all spots in the microarray.}

  Sources: Our own "research" and [2].
}

\section{Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{anonymize} \tab  -\cr
 \tab \code{as.RawData} \tab  -\cr
 \tab \code{getArea} \tab  -\cr
 \tab \code{getBackground} \tab  -\cr
 \tab \code{getBackgroundSD} \tab  -\cr
 \tab \code{getForeground} \tab  -\cr
 \tab \code{getForegroundSD} \tab  -\cr
 \tab \code{getForegroundSNR} \tab  -\cr
 \tab \code{\link[aroma:getRawData.QuantArrayData]{getRawData}} \tab Gets the raw intensites from the QuantArray data structure.\cr
 \tab \code{\link[aroma:getSpotPosition.QuantArrayData]{getSpotPosition}} \tab Gets physical positions of the spots.\cr
 \tab \code{\link[aroma:plotSpatial.QuantArrayData]{plotSpatial}} \tab Creates a spatial plot of a slide.\cr
 \tab \code{plotSpatial3d} \tab  -\cr
 \tab \code{\link[aroma:read.QuantArrayData]{read}} \tab Reads one or several QuantArray files into one QuantArrayData object.\cr
 \tab \code{setLayout} \tab  -\cr
 \tab \code{\link[aroma:write.QuantArrayData]{write}} \tab Write a QuantArray result data file.\cr
}


 \bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] QuantArray software by PerkinElmer Life Sciences,
      \url{http://lifesciences.perkinelmer.com/}.\cr
  [2] Packard BioScience,
      QuantArray Microarray Analysis Software Manual, 2001.\cr
}

\examples{
 \dontrun{
  # At the moment there is no QuantArray sample files in the package...
  # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
  # ...will use a GenePix sample file instead.
  qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  # Get the raw data; Rfg, Gfg, Rbg and Gbg
  raw <- getRawData(qa)

  # The the background corrected data
  ma <- getSignal(raw)

  # Plot M vs A with a lowess line through the data points
  subplots(4)
  plot(ma)
  plotSpatial(ma)
  boxplot(ma, groupBy="printtip")
 }
}

\note{
  The laser beam with wavelength 635nm is the red laser and excites the
  Cy5 dye, and the one with wavelength 532nm is the green laser which
  excites the Cy3 dye.
}
\keyword{classes}
