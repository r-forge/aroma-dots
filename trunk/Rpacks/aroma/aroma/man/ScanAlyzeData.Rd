%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyzeData.R
% 
% on Tue Jan 15 18:36:51 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ScanAlyzeData}
\docType{class}
\alias{ScanAlyzeData}


\title{The ScanAlyzeData class}

\usage{ScanAlyzeData(layout=NULL)}

\description{
 Package:  aroma \cr
\bold{Class ScanAlyzeData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{ScanAlyzeData}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{ScanAlyze20Data}\cr

public static class \bold{ScanAlyzeData}\cr
extends \link[aroma]{MicroarrayData}\cr



 Creates an empty ScanAlyzeData object.
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
 \tab \code{getArea} \tab  -\cr
 \tab \code{getBackground} \tab  -\cr
 \tab \code{getBgArea} \tab  -\cr
 \tab \code{getCircularity} \tab  -\cr
 \tab \code{getDiameter} \tab  -\cr
 \tab \code{getForeground} \tab  -\cr
 \tab \code{getGrid} \tab  -\cr
 \tab \code{\link[aroma:getRawData.ScanAlyzeData]{getRawData}} \tab Gets the raw (foreground and background) intensites.\cr
 \tab \code{getSpotColumn} \tab  -\cr
 \tab \code{\link[aroma:getSpotPosition.ScanAlyzeData]{getSpotPosition}} \tab Gets physical positions of the spots.\cr
 \tab \code{getSpotRow} \tab  -\cr
 \tab \code{normalizeGenewise} \tab  -\cr
 \tab \code{\link[aroma:plotSpatial.ScanAlyzeData]{plotSpatial}} \tab Creates a spatial plot of a slide.\cr
 \tab \code{plotSpatial3d} \tab  -\cr
 \tab \code{\link[aroma:read.ScanAlyzeData]{read}} \tab Reads one or several ScanAlyze files into a ScanAlyzeData object.\cr
 \tab \code{\link[aroma:write.ScanAlyzeData]{write}} \tab Write a ScanAlyze Results Data file.\cr
}


 \bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\section{Details}{
 A ScanAlyze file contains spot information for each spot on a single
 microarray slide. It consists of a header followed by a unspecified
 number of rows. The header contains a number of field labels. Each row
 starting with \code{SPOT} corresponds to one spot. Other rows starts
 with \code{HEADER} or \code{REMARK}. The spot fields are:\cr

  \item{SPOT}{Unique index of spot in file. "Counting starts with
    grid 1, moves along row 1 from column 1 until the last column, then
    advances to the next row; after all rows in grid 1 are assigned an
    index, counting proceeds to grid 2, etc.". Comment: This is the same
    way as GenePix and Spot indices the spots.}
  \item{GRID}{Grid number.}
  \item{TOP, LEFT, BOT, RIGHT}{Coordinates of the box containing spot
    ellipse, in image coordinates.}
  \item{ROW}{Row within grid.}
  \item{COL}{Column within grid.}
  \item{CH1I, CH2I}{Channel 1,2 - uncorrected mean of the foreground pixel
    intensites.}
  \item{CH1B, CH2B}{Channel 1,2 - median of the background pixel intensites.}
  \item{CH1BA, CH2BA}{Channel 1,2 - mean of the background pixel intensites.}
  \item{SPIX}{Number of pixels in the foreground, i.e. in the spot.}
  \item{BGPIX}{Number of pixels in the background.}
  \item{EDGE}{?}
  \item{RAT2}{Median of (CH2PI-CH2B)/(CH1PI-CH1B) where CH1PI and CH2PI are
    the single pixel intensities.}
  \item{MRAT}{?}
  \item{REGR}{?}
  \item{CORR}{The correlation between channel1 and channel2 pixels within
    the spot.}
  \item{LFRAT}{?}
  \item{CH1GTB1, CH2GTB1}{Fraction of pixels in the spot greater than the
    background (CH1B or CH2B).}
  \item{CH1GTB2, CH2GTB2}{Fraction of pixels in the spot greater than 1.5 times the background (CH1B or CH2B).}
  \item{CH1EDGEA, CH2EDGEA}{mean magnitude of the horizontal and vertical Sobel edge vectors within spot 1 and spot 2, respectively}
  \item{FLAG}{User defined spot flag (default 0).}
  \item{CH1KSD, CH2KSD}{The value of the Komogorov-Smirnov statistic [3] that assesses the likelihood that the spot pixel intensity distribution is drawn from the background distribution.}
  \item{CH1KSP, CH2KSP}{The actual probabilities of the above statistic.}
}


\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
 [1] ScanAlyze Software, Eisen Lab, Lawrence Berkeley National Lab, 2001,
  \url{http://rana.lbl.gov/EisenSoftware.htm}.\cr
 [2] ScanAlyze User Manual, Michael Eisen, Stanford University, 1999,
  \url{http://rana.lbl.gov/manuals/ScanAlyzeDoc.pdf}.\cr
 [3] 14.3 Are Two Distributions Different? a sample chapter from
     Numerical Recipes in C: The Art of Scientific Computing,
     1992, Cambridge University Press
  \url{http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf/c14-3.pdf}.\cr
}

\examples{
  sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))

  # Get the foreground and the background (and the layout)
  raw <- getRawData(sa)

  # The the background corrected data
  ma <- getSignal(raw, bgSubtract=FALSE)

  # Plot M vs A with a lowess line through the data points
  plot(ma, slide=1)
  lowessCurve(ma, lwd=2, gridwise=TRUE)
}

\note{
  The laser beam with wavelength 635nm is the red laser, and the one with
  wavelength 532nm is the green laser.
}
\keyword{classes}
