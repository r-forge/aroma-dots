%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenePixData}
\docType{class}
\alias{GenePixData}


 \title{The GenePixData class}

 \usage{GenePixData(layout=NULL, ...)}

 \description{
  Package:  aroma \cr
\bold{Class GenePixData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{GenePixData}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GenePixData}\cr
extends \link[aroma]{MicroarrayData}\cr


 }

 \arguments{
   \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
    slides in this data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{anonymize} \tab  -\cr
 \tab \code{append} \tab  -\cr
 \tab \code{as.RawData} \tab  -\cr
 \tab \code{getAbscent} \tab  -\cr
 \tab \code{\link[aroma:getArea.GenePixData]{getArea}} \tab Gets the area of the spot.\cr
 \tab \code{getArrayAspectRatio} \tab  -\cr
 \tab \code{getArrayBottom} \tab  -\cr
 \tab \code{getArrayHeight} \tab  -\cr
 \tab \code{getArrayLeft} \tab  -\cr
 \tab \code{getArrayRight} \tab  -\cr
 \tab \code{getArrayTop} \tab  -\cr
 \tab \code{getArrayWidth} \tab  -\cr
 \tab \code{getBackground} \tab  -\cr
 \tab \code{getBackgroundArea} \tab  -\cr
 \tab \code{\link[aroma:getBackgroundSD.GenePixData]{getBackgroundSD}} \tab Gets the standard deviation of the background pixels.\cr
 \tab \code{getBad} \tab  -\cr
 \tab \code{getBgArea} \tab  -\cr
 \tab \code{getCircularity} \tab  -\cr
 \tab \code{getFocusPosition} \tab  -\cr
 \tab \code{getForeground} \tab  -\cr
 \tab \code{\link[aroma:getForegroundSD.GenePixData]{getForegroundSD}} \tab Gets the standard deviation of the foreground pixels.\cr
 \tab \code{\link[aroma:getForegroundSE.GenePixData]{getForegroundSE}} \tab Gets the standard error of the foreground pixels.\cr
 \tab \code{getHeaderField} \tab  -\cr
 \tab \code{getHeaderFields} \tab  -\cr
 \tab \code{getLaserOnTime} \tab  -\cr
 \tab \code{getLaserPower} \tab  -\cr
 \tab \code{getLinesAveraged} \tab  -\cr
 \tab \code{getNotFound} \tab  -\cr
 \tab \code{getPixelSize} \tab  -\cr
 \tab \code{getPMTGain} \tab  -\cr
 \tab \code{\link[aroma:getRawData.GenePixData]{getRawData}} \tab Gets the raw intensites from the GPR structure.\cr
 \tab \code{getRgn R2} \tab  -\cr
 \tab \code{\link[aroma:getSpotPosition.GenePixData]{getSpotPosition}} \tab Gets physical positions of the spots.\cr
 \tab \code{getTemperature} \tab  -\cr
 \tab \code{getWavelengths} \tab  -\cr
 \tab \code{normalizeGenewise} \tab  -\cr
 \tab \code{\link[aroma:plotSpatial.GenePixData]{plotSpatial}} \tab Creates a spatial plot of a slide.\cr
 \tab \code{plotSpatial3d} \tab  -\cr
 \tab \code{\link[aroma:read.GenePixData]{read}} \tab Reads one or several GenePix files into one GenePixData object.\cr
 \tab \code{\link[aroma:readHeader.GenePixData]{readHeader}} \tab Reads a GenePix Results (GPR) file header.\cr
 \tab \code{setHeaderField} \tab  -\cr
 \tab \code{setLayout} \tab  -\cr
 \tab \code{\link[aroma:write.GenePixData]{write}} \tab Write a GenePix Results Data file.\cr
}


\bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

 }

 \details{
   For information about the fields in this class see Axon's
   specifications of the GenePix File Formats at
   \url{http://www.axon.com/GN_GenePix_File_Formats.html}.\cr
 }

 \section{Flagging feature-indicators}{
   GenePix Pro can \emph{flags} individual spots according to
   different criterias. The different flag statuses are
   \emph{Good} (has value 100 in the GPR structure),
   \emph{Bad} (-100), \emph{absent} (-75), and
   \emph{Not Found} (-50). Unflagged spots have value 0.
   A spot can only have one of these flags set at each time.

   The \emph{Absent} flag is set \emph{automatically} if
   there is a blank/empty entry in the GAL file (commonly
   reflected in the \code{\link{Layout}} object associated with this
   \code{\link{GenePixData}} object).

   The \emph{Not Found} flag is set \emph{automatically} if
   the spot/block alignment failed to find the spot. The
   position and diameter reported for such a spot are taken
   from the initial values obtained from the rough gridding.
   Note that the spot signals are still estimated.
   It may happen that the spot is dislocated, but that the
   user can visually pick it up and recenter the spot position
   and ask GenePix Pro to resegment the spot. The result
   may then be that GenePix flags the spot as \emph{Good}
   or \emph{Not Found} if it still finds that the spot is too
   weak for segmentation.

   The \emph{Bad} and \emph{Good} flags are set
   \emph{manually} by the user.
   A typical scenario is that GenePix identified the spots
   to have strong signals and segmented them well, but when
   the user looks at the image she or he sees that two spots
   are overlapping and have probably mixed their contents.
   Then the user flags that spot as \emph{Bad}.
   It can also be that there is a spot of interest, say a
   negative control, and GenePix flags it as \emph{Not Found},
   but the users looks at it and concludes that it is indeed
   a high quality negative control spot and therefore flags
   it as \emph{Good}.

   Spots that have been flagged \emph{Absent}, \emph{Bad},
   and \emph{Good} (not sure about this last one), will
   \emph{not} be relabel by GenePix. Thus, GenePix will
   only flag or unflag spots that are \emph{Not Found} or
   unflagged.

   Note that any flag may be set or unset manually after the
   automatic flagging has been done.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   GenePix File Formats,
   \url{http://www.axon.com/GN_GenePix_File_Formats.html} and
   \url{http://meetings.cshl.org/tgac/tgac/microarray.html}
 }

 \examples{
   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

   # Get the foreground and the background (and the layout)
   raw <- getRawData(gpr)

   # The the background corrected data
   ma <- getSignal(raw, bgSubtract=FALSE)

   # Plot M vs A with a lowess line through the data points
   plot(ma, slide=1)
   lowessCurve(ma, lwd=2, gridwise=TRUE)
 }

 \note{
   The laser beam with wavelength 635nm is the red laser and excites the
   Cy5 dye, and the one with wavelength 532nm is the green laser which
   excites the Cy3 dye.
 }
\keyword{classes}
