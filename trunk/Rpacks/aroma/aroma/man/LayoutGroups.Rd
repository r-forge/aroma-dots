%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  LayoutGroups.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{LayoutGroups}
\docType{class}
\alias{LayoutGroups}


\title{The LayoutGroups class}

\description{
 Package:  aroma \cr
\bold{Class LayoutGroups}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{LayoutGroups}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma]{GeneGroups}}, \emph{\link[aroma]{PlateGroups}}, \emph{\link[aroma]{PrintdipGroups}}, \emph{\link[aroma]{PrinttipGroups}}, \emph{\link[aroma]{SlideColumnGroups}}, \emph{\link[aroma]{SlideRowGroups}}, \emph{\link[aroma]{SuperGroups}}\cr

public abstract static class \bold{LayoutGroups}\cr
extends \link[R.oo]{Object}\cr



 The LayoutGroups class is a companion class for the Layout class.
}

\usage{LayoutGroups(layout=NULL)}

\arguments{
 \item{layout}{A \code{\link{Layout}} object specifying the layout of a set of
  microarrays.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{apply} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{getGroupValues} \tab  -\cr
 \tab \code{getLayout} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getSpots} \tab  -\cr
 \tab \code{getSpotValues} \tab  -\cr
 \tab \code{nbrOfGroups} \tab  -\cr
 \tab \code{seq} \tab  -\cr
 \tab \code{setLayout} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

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
  For more information see \code{\link{LayoutGroups}}.
}
\keyword{classes}
