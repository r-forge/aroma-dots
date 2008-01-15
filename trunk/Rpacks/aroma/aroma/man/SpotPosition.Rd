%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotPosition.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SpotPosition}
\docType{class}
\alias{SpotPosition}


\title{Class representing spot positions}

\description{
 Package:  aroma \cr
\bold{Class SpotPosition}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{SpotPosition}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SpotPosition}\cr
extends \link[R.oo]{Object}\cr



  A SpotPosition object contains physical (spatial) coordinates (in units
  of pixels) of the spots on a microarray slide.
}

\usage{SpotPosition(x=NULL, y=NULL)}

\arguments{
  \item{x}{An NxM \code{\link[base]{matrix}} with x positions, where N is the number of
    spots on each slide and M is the number of slides. Alternatively,
    another SpotPosition object can be given and then argument \code{y}
    is ignored.}
  \item{y}{An NxM \code{\link[base]{matrix}} with y positions of the same dimension as
    argument \code{x}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:as.character.SpotPosition]{as.character}} \tab Gets a character represenation of the spot positions.\cr
 \tab \code{\link[aroma:as.data.frame.SpotPosition]{as.data.frame}} \tab Gets a data frame represenation of the spot positions.\cr
 \tab \code{equals} \tab  -\cr
 \tab \code{extract} \tab  -\cr
 \tab \code{\link[aroma:getAspectRatio.SpotPosition]{getAspectRatio}} \tab Gets the aspect ratio of the spots.\cr
 \tab \code{\link[aroma:getBottomEdge.SpotPosition]{getBottomEdge}} \tab Gets the y-coordinate of the spot at the very bottom.\cr
 \tab \code{\link[aroma:getDistancesTo.SpotPosition]{getDistancesTo}} \tab Gets the Euclidean distance from one spot to a set of other spots.\cr
 \tab \code{\link[aroma:getLeftEdge.SpotPosition]{getLeftEdge}} \tab Gets the x-coordinate of the left most spot.\cr
 \tab \code{\link[aroma:getMaxHeight.SpotPosition]{getMaxHeight}} \tab Gets the maximum y-distance between two spots.\cr
 \tab \code{\link[aroma:getMaxWidth.SpotPosition]{getMaxWidth}} \tab Gets the maximum x-distance between two spots.\cr
 \tab \code{\link[aroma:getRightEdge.SpotPosition]{getRightEdge}} \tab Gets the x-coordinate of the right most spot.\cr
 \tab \code{\link[aroma:getTopEdge.SpotPosition]{getTopEdge}} \tab Gets the y-coordinate of the spot at the very top.\cr
 \tab \code{\link[aroma:nbrOfSlides.SpotPosition]{nbrOfSlides}} \tab Gets the number of slides.\cr
 \tab \code{\link[aroma:nbrOfSpots.SpotPosition]{nbrOfSpots}} \tab Gets the number of spots in each of the slides.\cr
 \tab \code{\link[aroma:plot.SpotPosition]{plot}} \tab Plots the spot positions in an xy scatter plot.\cr
 \tab \code{\link[aroma:points.SpotPosition]{points}} \tab Adds spot positions to an existing plot.\cr
 \tab \code{\link[aroma:size.SpotPosition]{size}} \tab Gets the number of spots in all slides together.\cr
 \tab \code{\link[aroma:write.SpotPosition]{write}} \tab Write a SpotPosition object to file.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
