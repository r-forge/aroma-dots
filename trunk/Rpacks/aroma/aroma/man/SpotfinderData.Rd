%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotfinderData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SpotfinderData}
\docType{class}
\alias{SpotfinderData}


\title{The SpotfinderData class}

\usage{SpotfinderData(layout=NULL)}

\description{
 Package:  aroma \cr
\bold{Class SpotfinderData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{SpotfinderData}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SpotfinderData}\cr
extends \link[aroma]{MicroarrayData}\cr



 Creates an SpotfinderData object. If the data frame \code{data} is empty or
 \code{\link[base]{NULL}}, the object will be empty.
}

\arguments{
  \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
   slides in this data set.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{append} \tab  -\cr
 \tab \code{as.RawData} \tab  -\cr
 \tab \code{extractLayout} \tab  -\cr
 \tab \code{getArea} \tab  -\cr
 \tab \code{getBackground} \tab  -\cr
 \tab \code{getBgArea} \tab  -\cr
 \tab \code{getCircularity} \tab  -\cr
 \tab \code{getForeground} \tab  -\cr
 \tab \code{getPerimeter} \tab  -\cr
 \tab \code{\link[aroma:getRawData.SpotfinderData]{getRawData}} \tab Gets the raw intensites from the SpotfinderData structure.\cr
 \tab \code{getSNR} \tab  -\cr
 \tab \code{normalizeGenewise} \tab  -\cr
 \tab \code{\link[aroma:plotSpatial.SpotfinderData]{plotSpatial}} \tab Creates a spatial plot of a slide.\cr
 \tab \code{plotSpatial3d} \tab  -\cr
 \tab \code{\link[aroma:read.SpotfinderData]{read}} \tab Reads several Spotfinder/TAV files into a SpotfinderData object.\cr
 \tab \code{\link[aroma:write.SpotfinderData]{write}} \tab Write a SpotfinderData object to file.\cr
}


\bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}

\section{Details}{
 A Spotfinder/TAV (TIGR Array Viewer) file is a headerless tab-delimited
 file with file extension \code{*.tav}, which the following contains 8 or
 17 columns:\cr

 The first eigth column do always exist:\cr
  \item{1}{spot row number. Range [0,ROWS] in Z+.}
  \item{2}{spot column number. Range [0,COLS] in Z+.}
  \item{3}{spot metarow/subrow number. Range [0,MROWS] in Z+.}
  \item{4}{spot metacolumn/subcolumn number. Range [0,MCOLS] in Z+.}
  \item{5}{spot row number in a block/grid. Range [0,GROWS] in Z+.}
  \item{6}{spot column number in a block/grid. Range [0,GCOLS] in Z+.}
  \item{7}{spot (integrated) intensity in channel A corrected for
    (median local) background. Note that saturated pixels are not used for calculations
    of intensities as they can skew the results. Range [0,MAXAREA*65535]
    in Z+.}
  \item{8}{spot (integrated) intensity in channel B corrected for
    (median local) background.  Note that saturated pixels are not used for calculations
    of intensities as they can skew the results. Range [0,MAXAREA*65535]
    in Z+.}

 Column 9-17 exists only in full TAV files:\cr
  \item{9}{spot mean ratio. Range [0,?] in R+. == (Col 7)/(Col 8), i.e. \bold{Redundant}.}

  \item{10}{spot total area in pixels. Range [0,MAXAREA] in Z+.}
  \item{11}{spot non-saturation factor. This measure shows the percentage
    of non-saturated pixels in the spot used for integration, i.e. values
    close to zero means a lot of saturation and values close to one means
    little saturation. Saturated pixels are not used for calculations of
    intensities as they can skew the results. Range [0,1] in R.}
  \item{12}{spot median ratio = (median(signalA)-median(bkgA))/(median(signalB)-median(bkgB)). Range [0,MAXRATIO] in R.}
  \item{13}{spot mode ratio = mode(signalA)-mode(nkgA))/(mode(signalB)-mode(bkgB)). Range [0,MAXRATIO] in R.}
  \item{14}{total spot background in channel A =(median(BkgA)* spotarea). This background was subtracted when the spot intensity was calculated. The background is calculated by measuring the average local background and multiplying it by the spot area. Range [0,MAXAREA*65535] in Z+.}
  \item{15}{total spot background in channel B =(median(BkgB)* spotarea). This background was subtracted when the spot intensity was calculated. The background is calculated by measuring the average local background and multiplying it by the spot area. Range [0,MAXAREA*65535] in Z+.}
  \item{16}{spot flag in channel A. This flag is set by QC filter. Range: \{A,B,C,X,Y,Z\}.}
  \item{17}{spot flag in channel B. This flag is set by QC filter.Range: \{A,B,C,X,Y,Z\}.}

  The foreground estimates for channel A (B) can be obtained as the sum of column 7 and 14 (column 8 and 14).
}

\section{Information about spot QC flags}{
  Flag values are generated based on next conditions:
   \item{A}{number of non-saturated pixels in spot is 0.
            \bold{Redundant}}
   \item{B}{number of non-saturated pixels in spot is between 0 and 50.
            \bold{Redundant}}
   \item{C}{number of non-saturated pixels in spot is more then 50.
            \bold{Redundant}}
   \item{X}{spot was detected and rejected based on spot shape and spot
            intensity relative to surrounding background.}
   \item{Y}{spot background is higher than spot intensity.
            \bold{Redundant}}
   \item{Z}{spot was not detected by the program.
            \bold{Redundant}}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  TIGR Spotfinder Software, The Institute for Genomic Research, USA,
  \url{http://www.tigr.org/software/tm4/}
}

\examples{\dontrun{
  sf <- SpotfinderData$read("sf111.tav", path=system.file("data-ex", package="aroma"))

  # Get the foreground and the background (and the layout)
  raw <- getRawData(spot)

  # The the background corrected data
  ma <- getSignal(raw, bgSubtract=FALSE)

  subplots(4, ncol=2)

  # Plot R vs G with a lowess line through the data points
  rg <- getRGData(ma)
  plot(rg)
  lowessCurve(rg, lwd=2, gridwise=TRUE)

  # Plot M vs A with a lowess line through the data points
  plot(ma)
  lowessCurve(ma, lwd=2, gridwise=TRUE)

  # Plot spatial
  plotSpatial(ma)
}}
\keyword{classes}
