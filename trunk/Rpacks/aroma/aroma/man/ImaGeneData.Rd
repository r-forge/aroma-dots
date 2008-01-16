%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ImaGeneData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ImaGeneData}
\docType{class}
\alias{ImaGeneData}


\title{The ImaGeneData class}

\usage{ImaGeneData(layout=NULL)}

\description{
 Package:  aroma \cr
\bold{Class ImaGeneData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{MicroarrayData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{ImaGeneData}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ImaGeneData}\cr
extends \link[aroma]{MicroarrayData}\cr



 Creates an empty ImaGeneData object.
}

\arguments{
  \item{layout}{A \code{\link{Layout}} object specifying the spot layout of the
   slides in this data set.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{anonymize} \tab  -\cr
 \tab \code{as.RawData} \tab  -\cr
 \tab \code{getAbscent} \tab  -\cr
 \tab \code{getArea} \tab  -\cr
 \tab \code{getBgArea} \tab  -\cr
 \tab \code{getCircularity} \tab  -\cr
 \tab \code{\link[aroma:getRawData.ImaGeneData]{getRawData}} \tab Gets the raw intensites from two ImaGeneData objects.\cr
 \tab \code{\link[aroma:getSpotPosition.ImaGeneData]{getSpotPosition}} \tab Gets physical positions of the spots.\cr
 \tab \code{plotSpatial} \tab  -\cr
 \tab \code{plotSpatial3d} \tab  -\cr
 \tab \code{\link[aroma:read.ImaGeneData]{read}} \tab Reads an ImaGene result file.\cr
 \tab \code{setLayout} \tab  -\cr
 \tab \code{\link[aroma:write.ImaGeneData]{write}} \tab Write an ImaGene result data file.\cr
}


 \bold{Methods inherited from MicroarrayData}:\cr
addFlag, append, applyGenewise, applyGroupwise, applyPlatewise, applyPrintdipwise, applyPrinttipwise, as.character, as.data.frame, boxplot, clearCache, clearFlag, createColors, dataFrameToList, equals, extract, getBlank, getCache, getChannelNames, getColors, getExcludedSpots, getExtra, getExtreme, getFieldNames, getFlag, getInclude, getLabel, getLayout, getProbeWeights, getSignalWeights, getSlideNames, getSlidePairs, getSpotPosition, getSpotValue, getTreatments, getView, getWeights, getWeightsAsString, hasExcludedSpots, hasLayout, hasProbeWeights, hasSignalWeights, hasWeights, highlight, hist, isFieldColorable, keepSlides, keepSpots, listFlags, lowessCurve, nbrOfDataPoints, nbrOfFields, nbrOfSlides, nbrOfSpots, nbrOfTreatments, normalizePlatewise, normalizePrintorder, normalizeQuantile, plot, plotDensity, plotGene, plotPrintorder, plotReplicates, plotSpatial, plotSpatial3d, plotXY, points, putGene, putSlide, qqnorm, quantile, range, range2, read, readHeader, readToList, removeSlides, removeSpots, resetProbeWeights, resetSignalWeights, select, seq, setCache, setExcludedSpots, setExtra, setFlag, setLabel, setLayout, setProbeWeights, setSignalWeights, setSlideNames, setTreatments, setView, setWeights, size, str, subplots, summary, text, updateHeader, validateArgumentChannel, validateArgumentChannels, validateArgumentGroupBy, validateArgumentSlide, validateArgumentSlides, validateArgumentSpotIndex, validateArgumentWeights, write, writeHeader

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  BioDiscovery's ImaGene software,
  \url{http://www.biodiscovery.com/imagene.asp}
}

\examples{
######################################################################
# Example that reads two ImaGene data files, which are single-channel
# and represent one slide, and combines them into a RawData object.
######################################################################
igG <- ImaGeneData$read("imagene234.cy3", path=system.file("data-ex", package="aroma"))
igR <- ImaGeneData$read("imagene234.cy5", path=system.file("data-ex", package="aroma"))

# Get the raw data; Rfg, Gfg, Rbg and Gbg
raw <- getRawData(igR, igG)

# The the background corrected data
ma <- getSignal(raw)

# Plot M vs A with a lowess line through the data points
subplots(4)
plot(ma)
plotSpatial(ma)
boxplot(ma, groupBy="printtip")


######################################################################
# Example that reads multiple ImaGene data files, which all are
# single-channel and represent several slide, and then combines them
# into a RawData object of several slides.
######################################################################

# Function that reads several (single-channel) ImaGene data files
# and combines them into a multi-slide two-channel RawData object
myImaGeneReadAll <- function(slides, extensions=c("cy3","cy5"), ...) {
  # An empty RawData object to store all slides
  raw <- RawData()
  for (slide in slides) {
    igCy3 <- ImaGeneData$read(paste(slide, extensions[1], sep="."), ...)
    igCy5 <- ImaGeneData$read(paste(slide, extensions[2], sep="."), ...)
    tmp <- getRawData(igCy3, igCy5)
    append(raw, tmp)
  }

  # Return the RawData object containing all slides
  raw;
} # myImaGeneReadAll()


# "Faking" three slides...
slides <- rep("imagene234", 3)

# ...and reads them.
raw <- myImaGeneReadAll(slides, path=system.file("data-ex", package="aroma"))
}

\note{
  The laser beam with wavelength 635nm is the red laser and excites the
  Cy5 dye, and the one with wavelength 532nm is the green laser which
  excites the Cy3 dye.
}
\keyword{classes}
