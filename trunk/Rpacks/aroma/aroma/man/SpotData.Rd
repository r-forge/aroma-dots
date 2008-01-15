%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SpotData}
\docType{class}
\alias{SpotData}


 \title{The SpotData class}

 \usage{SpotData(layout=NULL)}

 \description{
  Package:  aroma \cr
\bold{Class SpotData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{SpotData}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SpotData}\cr
extends \link[aroma]{MicroarrayData}\cr



  Creates an SpotData object. If the data frame \code{data} is empty or
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
 \tab \code{getArrayAspectRatio} \tab  -\cr
 \tab \code{getArrayBottom} \tab  -\cr
 \tab \code{getArrayHeight} \tab  -\cr
 \tab \code{getArrayLeft} \tab  -\cr
 \tab \code{getArrayRight} \tab  -\cr
 \tab \code{getArrayTop} \tab  -\cr
 \tab \code{getArrayWidth} \tab  -\cr
 \tab \code{getBackground} \tab  -\cr
 \tab \code{getBgArea} \tab  -\cr
 \tab \code{getCircularity} \tab  -\cr
 \tab \code{getForeground} \tab  -\cr
 \tab \code{\link[aroma:getForegroundSD.SpotData]{getForegroundSD}} \tab Gets (an approximation of) the standard deviation of the foreground pixels.\cr
 \tab \code{\link[aroma:getForegroundSE.SpotData]{getForegroundSE}} \tab Gets the standard error of the foreground pixels.\cr
 \tab \code{getPerimeter} \tab  -\cr
 \tab \code{\link[aroma:getRawData.SpotData]{getRawData}} \tab Gets the raw intensites from the SpotData structure.\cr
 \tab \code{getSNR} \tab  -\cr
 \tab \code{\link[aroma:getSpotPosition.SpotData]{getSpotPosition}} \tab Gets physical positions of the spots.\cr
 \tab \code{log} \tab  -\cr
 \tab \code{normalizeGenewise} \tab  -\cr
 \tab \code{\link[aroma:plotSpatial.SpotData]{plotSpatial}} \tab Creates a spatial plot of a slide.\cr
 \tab \code{plotSpatial3d} \tab  -\cr
 \tab \code{\link[aroma:read.SpotData]{read}} \tab Reads several Spot files into a SpotData object.\cr
 \tab \code{readOneFile} \tab  -\cr
 \tab \code{\link[aroma:write.SpotData]{write}} \tab Write a SpotData object to file.\cr
}


\bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

 }

 \section{Details}{
  A Spot file contains spot information for each spot on a single microarray slide. It consists of a header followed by a unspecified number of rows. The header contains 1+30 labels, and each row contains 31 fields. Each row corresponds to one spot. The fields are:\cr

   \item{<NO NAME>}{row number}

   \item{indexs}{spot number on slide. Range [0,N] in N.}

   \item{grid.r}{grid row number. Range [1,GR] in Z+.}
   \item{grid.c}{grid column number. Range [1,GC] in Z+.}
   \item{spot.r}{spot row (within grid) number. Range [1,SR] in Z+.}
   \item{spot.c}{spot column (within grid) number. Range [1,SC] in Z+.}

   \item{area}{the number of foreground pixels. Range [0,MAXAREA] in N}

   \item{Gmean}{the average of the foreground pixel values. Range [0,65535] in R}
   \item{Gmedian}{the median of the foreground pixel values. Range [0,65535] in N}
   \item{GIQR}{the inter quartile range (a robust measure of spread) of the logged foregroud pixel values. Range [0,16]+Inf in {R,NA}}
   \item{Rmean}{the average of the foreground pixel values. Range [0,65535] in R}
   \item{Rmedian}{the median of the foreground pixel values. Range [0,65535] in N}
   \item{RIQR}{the inter quartile range (a robust measure of spread) of the logged foregroud pixel values. Range [0,16]+Inf in {R,NA}}

   \item{bgGmean}{mean green background intesity. Range [0,65535] in R}
   \item{bgGmed}{median green background intesity. Range [0,65535] in N}
   \item{bgGSD}{standard deviation for the green background. Range [0,65535]+Inf in R}
   \item{bgRmean}{mean red background intesity. Range [0,65535] in R}
   \item{bgRmed}{median red background intesity. Range [0,65535] in N}
   \item{bgRSD}{standard deviation for the red background. Range [0,65535]+Inf in R}

   \item{valleyG}{the background intesity estimate from the local background valley method S.valley. Range [0,65535] in N}
   \item{valleyR}{the background intesity estimate from the local background valley method S.valley. Range [0,65535] in N}

   \item{morphG}{green background estimate using morphological opening (erosion-dilation). Range [0,65535] in N}
   \item{morphG.erode}{green background estimate using morphological erosion. Range [0,65535] in N}
   \item{morphG.close.open}{green background estimate using morphological closing-opening (dilation-erosion-dilation). Range [0,65535] in N}
   \item{morphR}{red background estimate using morphological opening (erosion-dilation). Range [0,65535] in N}
   \item{morphR.erode}{red background estimate using morphological erosion. Range [0,65535] in N}
   \item{morphR.close.open}{red background estimate using morphological closing-opening (dilation-erosion-dilation). Range [0,65535] in N}

   \item{logratio}{== log((Rmedian-morphR)/(Gmedian-morphG), base=2), i.e. \bold{Redundant}.}
   \item{perimeter}{== 2*sqrt(pi*area/circularity), i.e. \bold{Redundant}.}
   \item{circularity}{Shape of spot defined as 4*pi*area/perimeter**2.}
   \item{badspot}{If the spot area is greater than product of the horizontal and the vertical average spot separations, equal to \code{1}, otherwise \code{0}.}
 }

 \section{About IQR}{
   The interquartile range (IQR) is the distance between the 75\%
   \code{\link[base]{quantile}} (percentile) and the 25\% quantile.
   In words, IQR is the range of the mid 50\%. Thus, no outliers are
   included in the measure, which is why we say it is a robust measure.
   For norammly distributed data \eqn{IQR = 1.35*\sigma}, where \eqn{\sigma}
   is the standard deviation.
 }

 \section{More on background estimates}{
   The Spot software provides several different kinds of background estimates where three of them are
   based on morphological methods. For all of these methods, the signal selected to be the background
   signal is the pixel value at the center of the spot \emph{after} applying the morphological transform
   using a square mask with side 2.5 times the average distance between two spots.
   The first and also the simpliest transform (\code{morph.erode}) performs a single \emph{erosion} step.
   The second transform (\code{morph}) performs an \emph{opening}, which is an \emph{erosion} followed by
   a \emph{dilution}.
   The third transform (\code{morph.close.open}) performs a \emph{closing} followed by an \emph{opening},
   which is the same as doing a \emph{dilution}, then an \emph{erosion} and a \emph{dilution} again.
   As the names of the steps indicate, an erosion makes the signal smaller and the dilution the signal
   larger. Hence, background estimated based on these three methods can always be ordered as
   \code{morph.erode <= morph <= morph.close.open}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   Spot Software package by CSIRO, Australia,
   \url{http://www.cmis.csiro.au/iap/spot.htm}

   Spot: Description of Output, 2003
   \url{http://www.cmis.csiro.au/iap/Spot/spotoutput.htm}

   Y.H. Yang, M. Buckley, S. Dudoit, T. Speed,
   Comparison of methods for image analysis on cDNA microarray data,
   Tech. Report 584, Nov 2000.
   \url{http://www.stat.berkeley.edu/users/terry/zarray/Html/image.html}
 }

 \examples{
   spot <- SpotData$read("spot123.spot", path=system.file("data-ex", package="aroma"))

   # Get the foreground and the background (and the layout)
   raw <- getRawData(spot)

   # The the background corrected data
   ma <- getSignal(raw, bgSubtract=FALSE)

   subplots(4, ncol=2)

   # Plot R vs G with a lowess line through the data points
   rg <- as.RGData(ma)
   plot(rg)
   lowessCurve(rg, lwd=2, gridwise=TRUE)

   # Plot M vs A with a lowess line through the data points
   plot(ma)
   lowessCurve(ma, lwd=2, gridwise=TRUE)

   # Plot spatial
   plotSpatial(ma)
 }
\keyword{classes}
