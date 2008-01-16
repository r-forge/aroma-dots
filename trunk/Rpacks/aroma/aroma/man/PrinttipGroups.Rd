%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PrinttipGroups.R
% 
% on Tue Jan 15 18:34:52 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{PrinttipGroups}
\docType{class}
\alias{PrinttipGroups}


\title{The PrinttipGroups class}

\description{
 Package:  aroma \cr
\bold{Class PrinttipGroups}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{LayoutGroups}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{PrinttipGroups}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{PrinttipGroups}\cr
extends \emph{\link[aroma]{LayoutGroups}}\cr


}

\usage{PrinttipGroups(layout=NULL, groups=NULL)}

\arguments{
 \item{layout}{A \code{\link{Layout}} object specifying the layout of a set of
  microarrays.}
 \item{groups}{A \code{\link[base]{list}} of length equal to the number of printtips and
  each element contains spot indices corresponding to each printtip. If
  \code{\link[base]{NULL}}, the printtips are obtained from the Layout object.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
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
  For more information see \code{\link{PrinttipGroups}}.
}
\keyword{classes}
