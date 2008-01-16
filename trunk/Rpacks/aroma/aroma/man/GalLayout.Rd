%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GalLayout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GalLayout}
\docType{class}
\alias{GalLayout}


\title{Class representing GenePix Array List (GAL) data}

\description{
 Package:  aroma \cr
\bold{Class GalLayout}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Layout}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{GalLayout}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GalLayout}\cr
extends \link[aroma]{Layout}\cr


}

\usage{GalLayout(...)}

\arguments{
  \item{...}{Arguments accepted by the Layout constructor.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.data.frame} \tab  -\cr
 \tab \code{getBlock} \tab  -\cr
 \tab \code{getColumn} \tab  -\cr
 \tab \code{getFieldNames} \tab  -\cr
 \tab \code{getGridColumn} \tab  -\cr
 \tab \code{getGridRow} \tab  -\cr
 \tab \code{getRow} \tab  -\cr
 \tab \code{nbrOfFields} \tab  -\cr
 \tab \code{\link[aroma:read.GalLayout]{read}} \tab Reads a GenePix Array List (GAL) file.\cr
 \tab \code{write} \tab  -\cr
}


 \bold{Methods inherited from Layout}:\cr
anonymize, as.character, as.data.frame, as, as.Layout, equals, fromDataFrame, getBlank, getGeneGroups, getId, getIndex, getIndices, getLayoutGroupsByName, getLocation, getName, getPlate, getPlateGroups, getPlateNumber, getPosition, getPrintdipGroups, getPrintorder, getPrinttipGroups, getSlideColumnGroups, getSlideRowGroups, gridSize, hasIds, hasNames, hasPlates, indexOf, nbrOfColumns, nbrOfGrids, nbrOfPlates, nbrOfReplicates, nbrOfRows, nbrOfSpots, put, read, rename, set, setGeneGroups, setId, setName, setPlate, setPlateGroups, setPrintdipGroups, setPrintorder, setPrinttipGroups, setSlideColumnGroups, setSlideRowGroups, size, toXYMatrix, write

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
 \dontrun{
  layout <- GalLayout$read("fish.gal")
  genes <- getGeneGroups(layout)
  print(genes)
 }

 \dontrun{
  # Gives:
   [1] "GeneGroups: 7769 groups defined. 7680 genes with 1 replicates, 75
   genes with 3 replicates, 5 genes with 6 replicates, 4 genes with 16
   replicates, 4 genes with 64 replicates, 1 genes (3XSSC) with 193
   replicates."
 }
}

\seealso{
 The superclass \code{\link{Layout}}.
}

\references{
 The GenePix Array List (GAL) file format,
 \url{http://www.axon.com/GN_GenePix_File_Formats.html}.
}
\keyword{classes}
