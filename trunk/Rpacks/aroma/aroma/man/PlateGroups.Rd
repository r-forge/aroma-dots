%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  LayoutGroups.R
% 
% on Tue Jan 15 18:33:05 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{PlateGroups}
\docType{class}
\alias{PlateGroups}


\title{The PlateGroups class}

\description{
 Package:  aroma \cr
\bold{Class PlateGroups}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{LayoutGroups}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{PlateGroups}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{PlateGroups}\cr
extends \emph{\link[aroma]{LayoutGroups}}\cr


}

\usage{PlateGroups(layout=NULL, plateDef=NULL)}

\arguments{
 \item{layout}{A \code{\link{Layout}} object specifying the layout of a set of
  microarrays.}
 \item{plateDef}{
  If a scalar (\code{\link[base]{numeric}} of length one), each plate group is assumed to
  contain that number of spots/clones.
  Can also be a \code{\link[base]{data.frame}}.
  If \code{\link[base]{NULL}}, the plates are obtained from the Layout
  object, which then is required to contain plate information.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{getFirst} \tab  -\cr
 \tab \code{getLast} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getPrintorderIndices} \tab  -\cr
 \tab \code{getSizes} \tab  -\cr
 \tab \code{getSpots} \tab  -\cr
 \tab \code{nbrOfGroups} \tab  -\cr
}


 \bold{Methods inherited from LayoutGroups}:\cr
apply, as.character, getGroupValues, getLayout, getNames, getSpots, getSpotValues, nbrOfGroups, seq, setLayout

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))

# Extract the printtip groups
printtips <- PrinttipGroups(layout)
print(printtips)

# Extract the plate groups (using getPlate(layout))
plates <- PlateGroups(layout)
print(plates)

# Extract the (within-slide) replicate groups (using getName/ID(layout))
replicates <- GeneGroups(layout)
print(replicates)

# Define the (within-slide) replicate groups
duplicates <- GeneGroups(layout, "neighboring-pairs")
print(duplicates)
}

\seealso{
  For more information see \code{\link{PlateGroups}}.
}
\keyword{classes}
