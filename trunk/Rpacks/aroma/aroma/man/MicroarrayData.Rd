%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{MicroarrayData}
\docType{class}
\alias{MicroarrayData}


\title{The MicroarrayData class}

\description{
 Package:  aroma \cr
\bold{Class MicroarrayData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{MicroarrayData}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{BMAData}, \link[aroma]{GenePixData}, \link[aroma]{ImaGeneData}, \link[aroma]{MAData}, \link[aroma]{QuantArrayData}, \link[aroma]{RawData}, \link[aroma]{RGData}, \link[aroma]{ScanAlyze20Data}, \link[aroma]{ScanAlyzeData}, \link[aroma]{SpotData}, \link[aroma]{SpotfinderData}, \link[aroma]{TMAData}\cr

public static class \bold{MicroarrayData}\cr
extends \link[R.oo]{Object}\cr



 This class is \emph{abstract} can not be instanciated.
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
   \tab \code{layout} \tab The layout structure of the microarray slide(s). \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:append.MicroarrayData]{append}} \tab Appends another MicroarrayData object.\cr
 \tab \code{applyGenewise} \tab  -\cr
 \tab \code{applyGroupwise} \tab  -\cr
 \tab \code{applyPlatewise} \tab  -\cr
 \tab \code{applyPrintdipwise} \tab  -\cr
 \tab \code{applyPrinttipwise} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma:as.data.frame.MicroarrayData]{as.data.frame}} \tab Converts the object to a data frame.\cr
 \tab \code{\link[aroma:boxplot.MicroarrayData]{boxplot}} \tab Plots a boxplot.\cr
 \tab \code{clearCache} \tab  -\cr
 \tab \code{\link[aroma:equals.MicroarrayData]{equals}} \tab Compares a MicroarrayData object with another Object.\cr
 \tab \code{\link[aroma:extract.MicroarrayData]{extract}} \tab Gets a subset of data specified by fields, slides and/or spot indices.\cr
 \tab \code{getBlank} \tab  -\cr
 \tab \code{getCache} \tab  -\cr
 \tab \code{\link[aroma:getChannelNames.MicroarrayData]{getChannelNames}} \tab Gets the names of the channels.\cr
 \tab \code{getColors} \tab  -\cr
 \tab \code{getExtreme} \tab  -\cr
 \tab \code{\link[aroma:getFieldNames.MicroarrayData]{getFieldNames}} \tab Gets the names of all fields containing spot specific data.\cr
 \tab \code{getInclude} \tab  -\cr
 \tab \code{\link[aroma:getLabel.MicroarrayData]{getLabel}} \tab Gets the label of one field.\cr
 \tab \code{\link[aroma:getLayout.MicroarrayData]{getLayout}} \tab Gets the layout.\cr
 \tab \code{getProbeWeights} \tab  -\cr
 \tab \code{getSignalWeights} \tab  -\cr
 \tab \code{\link[aroma:getSlideNames.MicroarrayData]{getSlideNames}} \tab Gets the names of the slides.\cr
 \tab \code{getSlidePairs} \tab  -\cr
 \tab \code{getSpotPosition} \tab  -\cr
 \tab \code{getSpotValue} \tab  -\cr
 \tab \code{getTreatments} \tab  -\cr
 \tab \code{getView} \tab  -\cr
 \tab \code{getWeights} \tab  -\cr
 \tab \code{getWeightsAsString} \tab  -\cr
 \tab \code{\link[aroma:hasLayout.MicroarrayData]{hasLayout}} \tab Checks if the layout has been specified.\cr
 \tab \code{hasProbeWeights} \tab  -\cr
 \tab \code{hasSignalWeights} \tab  -\cr
 \tab \code{hasWeights} \tab  -\cr
 \tab \code{\link[aroma:highlight.MicroarrayData]{highlight}} \tab Highlight spots in last plot.\cr
 \tab \code{\link[aroma:hist.MicroarrayData]{hist}} \tab Plots a histogram.\cr
 \tab \code{isFieldColorable} \tab  -\cr
 \tab \code{\link[aroma:keepSlides.MicroarrayData]{keepSlides}} \tab Remove all but the specified slides.\cr
 \tab \code{\link[aroma:keepSpots.MicroarrayData]{keepSpots}} \tab Remove all but the specified spots from all slides.\cr
 \tab \code{\link[aroma:lowessCurve.MicroarrayData]{lowessCurve}} \tab Draws one or more lowess curves through the data in last plot.\cr
 \tab \code{nbrOfDataPoints} \tab  -\cr
 \tab \code{\link[aroma:nbrOfFields.MicroarrayData]{nbrOfFields}} \tab Gets the number of fields.\cr
 \tab \code{\link[aroma:nbrOfSlides.MicroarrayData]{nbrOfSlides}} \tab Gets the number of slides.\cr
 \tab \code{\link[aroma:nbrOfSpots.MicroarrayData]{nbrOfSpots}} \tab Gets the number of spots in each of the slides.\cr
 \tab \code{nbrOfTreatments} \tab  -\cr
 \tab \code{\link[aroma:normalizePlatewise.MicroarrayData]{normalizePlatewise}} \tab Normalization performed plate by plate.\cr
 \tab \code{normalizePrintorder} \tab  -\cr
 \tab \code{normalizeQuantile} \tab  -\cr
 \tab \code{\link[aroma:plot.MicroarrayData]{plot}} \tab Plots spatial representation of a microarray spot statistics.\cr
 \tab \code{\link[aroma:plotDensity.MicroarrayData]{plotDensity}} \tab Plots the empirical density distribution of a field.\cr
 \tab \code{plotGene} \tab  -\cr
 \tab \code{\link[aroma:plotPrintorder.MicroarrayData]{plotPrintorder}} \tab Plots the data as a time series in the order it was printed.\cr
 \tab \code{plotReplicates} \tab  -\cr
 \tab \code{\link[aroma:plotSpatial.MicroarrayData]{plotSpatial}} \tab Plots a spatial representation of one field.\cr
 \tab \code{\link[aroma:plotSpatial3d.MicroarrayData]{plotSpatial3d}} \tab Plots a 3-dimensional spatial representation of a field.\cr
 \tab \code{\link[aroma:plotXY.MicroarrayData]{plotXY}} \tab Plots a scatter plot of two fields.\cr
 \tab \code{points} \tab  -\cr
 \tab \code{qqnorm} \tab  -\cr
 \tab \code{quantile} \tab  -\cr
 \tab \code{range} \tab  -\cr
 \tab \code{range2} \tab  -\cr
 \tab \code{\link[aroma:read.MicroarrayData]{read}} \tab Reads microarray data generated by another software package.\cr
 \tab \code{readHeader} \tab  -\cr
 \tab \code{\link[aroma:removeSlides.MicroarrayData]{removeSlides}} \tab Remove specified slides.\cr
 \tab \code{\link[aroma:removeSpots.MicroarrayData]{removeSpots}} \tab Remove specified spots.\cr
 \tab \code{resetProbeWeights} \tab  -\cr
 \tab \code{resetSignalWeights} \tab  -\cr
 \tab \code{\link[aroma:select.MicroarrayData]{select}} \tab Selects rows with certain criteria.\cr
 \tab \code{seq} \tab  -\cr
 \tab \code{setCache} \tab  -\cr
 \tab \code{\link[aroma:setLabel.MicroarrayData]{setLabel}} \tab Sets the label of one field.\cr
 \tab \code{\link[aroma:setLayout.MicroarrayData]{setLayout}} \tab Sets the layout.\cr
 \tab \code{\link[aroma:setProbeWeights.MicroarrayData]{setProbeWeights}} \tab Sets probe weights on one or several slides.\cr
 \tab \code{setSignalWeights} \tab  -\cr
 \tab \code{\link[aroma:setSlideNames.MicroarrayData]{setSlideNames}} \tab Sets the names of the slides.\cr
 \tab \code{setTreatments} \tab  -\cr
 \tab \code{setView} \tab  -\cr
 \tab \code{setWeights} \tab  -\cr
 \tab \code{\link[aroma:size.MicroarrayData]{size}} \tab Gets the number of spots in all slides together.\cr
 \tab \code{\link[aroma:str.MicroarrayData]{str}} \tab Compactly Display the Structure of a MicroarrayData object.\cr
 \tab \code{\link[aroma:subplots.MicroarrayData]{subplots}} \tab Prepare for a grid of subplots based on the number of slides.\cr
 \tab \code{\link[aroma:summary.MicroarrayData]{summary}} \tab Gets summary statistics for the fields.\cr
 \tab \code{\link[aroma:text.MicroarrayData]{text}} \tab Puts labels to the spots in last plot.\cr
 \tab \code{updateHeader} \tab  -\cr
 \tab \code{validateArgumentChannel} \tab  -\cr
 \tab \code{validateArgumentChannels} \tab  -\cr
 \tab \code{\link[aroma:validateArgumentGroupBy.MicroarrayData]{validateArgumentGroupBy}} \tab Validates the argument 'groupBy'.\cr
 \tab \code{\link[aroma:validateArgumentSlide.MicroarrayData]{validateArgumentSlide}} \tab Validates the argument 'slide'.\cr
 \tab \code{\link[aroma:validateArgumentSlides.MicroarrayData]{validateArgumentSlides}} \tab Validates the argument 'slides'.\cr
 \tab \code{validateArgumentSpotIndex} \tab  -\cr
 \tab \code{\link[aroma:validateArgumentWeights.MicroarrayData]{validateArgumentWeights}} \tab Validates the argument 'weights'.\cr
 \tab \code{\link[aroma:write.MicroarrayData]{write}} \tab Write a MicroarrayData object to file.\cr
 \tab \code{writeHeader} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
