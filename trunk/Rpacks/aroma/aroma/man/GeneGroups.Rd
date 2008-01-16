%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GeneGroups.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GeneGroups}
\docType{class}
\alias{GeneGroups}


\title{The GeneGroups class}

\description{
 Package:  aroma \cr
\bold{Class GeneGroups}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{LayoutGroups}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{GeneGroups}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{GeneGroups}\cr
extends \emph{\link[aroma]{LayoutGroups}}\cr


}

\usage{GeneGroups(layout=NULL, specification=NULL)}

\arguments{
  \item{layout}{A Layout object.}
  \item{specification}{(optional) \code{\link[base]{character}} string specifying what type of
    within-slide replication the slides are printed with. If \code{\link[base]{NULL}},
    the within-slide replication will be extracted from the names or the
    ids of the Layout object. Possible values are
    \code{"none"} (no within-slide replications),
    \code{"neighboring-pairs"} (duplicated horizontal neighboring pairs) and
    \code{"top-bottom"} (duplicated pairs over replicated tiles).
  }
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma:getDuplicates.GeneGroups]{getDuplicates}} \tab Gets all duplicated spots.\cr
 \tab \code{\link[aroma:getNames.GeneGroups]{getNames}} \tab Gets the names of the unique gene names.\cr
 \tab \code{\link[aroma:getQuadruplicates.GeneGroups]{getQuadruplicates}} \tab Gets all quadruplicated spots.\cr
 \tab \code{\link[aroma:getReplicates.GeneGroups]{getReplicates}} \tab Gets all spots that exist with a certain number of replicates.\cr
 \tab \code{\link[aroma:getSizes.GeneGroups]{getSizes}} \tab Gets the number of replicates for each gene.\cr
 \tab \code{\link[aroma:getSpots.GeneGroups]{getSpots}} \tab Gets a list of items each containing spot indices for a group.\cr
 \tab \code{\link[aroma:getTriplicates.GeneGroups]{getTriplicates}} \tab Gets all triplicated spots.\cr
 \tab \code{\link[aroma:nbrOfGroups.GeneGroups]{nbrOfGroups}} \tab Gets the number of unique genes.\cr
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
\keyword{classes}
