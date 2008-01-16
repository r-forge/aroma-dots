%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Layout}
\docType{class}
\alias{Layout}


\title{The Layout class}

\description{
 Package:  aroma \cr
\bold{Class Layout}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Layout}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{GalLayout}\cr

public static class \bold{Layout}\cr
extends \link[R.oo]{Object}\cr



 The Layout class describes the layout of a microarray slide,
 such as the number of spots and the layout of the grids etc.
 A microarray slide is layout in a number of grids, where each grid
 is refered to by its row and its column. Within each grid, each spot
 is refered to by its row and its column \emph{within that grid}. For
 example, a microarray slide with 4*3 grids and where each grid has
 12*10 spots, has in total 4*3*12*10 = 12*120 = 1440 spots. The spot
 in the top left corner is located at grid (1,1) and spot (1,1) and
 we say it has the \emph{location} (1,1,1,1). In the above example,
 the spot in the lower right corner of the slide has the location
 (4,3,12,10). Further more, each spot on a microarray slide has a
 unique \emph{index}. This index starts counting starts with grid 1,
 moves along row 1 from column 1 until
 the last column, then advances to the next row. After all rows in
 grid 1 are indexed, counting proceeds to grid 2 and so on. This way of
 indexing the spots is used by for instance GenePix, Spot and ScanAlyze.
 Continuing the example above, the spot at location (1,1,1,1) has index
 1, the spot at location (4,3,12,10) has index 1440. The spot at location
 (2,2,11,3) has index ((2-1)*3+(2-1))*12*10+(11-1)*10+3 = 483.
 To get the \emph{index} from a \emph{location} one do
 \code{getIndex(layout, 2,2,11,3)} where
 \code{layout <- Layout(4,3,12,10)}. To get the \emph{location} from a
 \emph{index} one do \code{getLocation(layout, 483)}, which
 gives the vector (2, 2, 11, 3). To get the \emph{position}, i.e. the
 overall row and column of a spot on the microarray slide one can do
 \code{getPosition(483)} which give the position (23,13).\cr

 All spots are refered to by their unique \emph{indices}.
}

\usage{Layout(ngrid.r=0, ngrid.c=0, nspot.r=0, nspot.c=0, geneSpotMap=NULL, name=NULL, id=NULL, plate=NULL, printorder=NULL)}

\arguments{
  \item{nspot.r}{The number of rows of spots per grid.
        This first argument can also be a list containing
        the fields nspot.r, nspot.c, ngrid.r, and ngrid.c.}
  \item{nspot.c}{The number of columns of spots per grid.}
  \item{ngrid.r}{The number of rows of grids per slide.}
  \item{ngrid.c}{The number of columns of grids per slide.}
  \item{name}{A \code{\link[base]{vector}} of names for the spot. \bold{optional}.}
  \item{id}{A \code{\link[base]{vector}} of ids for the spot. \bold{optional}.}
  \item{printorder}{A \code{\link[base]{matrix}} or a \code{\link[base]{character}} string specifying the
        order that spots where printed. For more details
        \code{\link[aroma:setPrintorder.Layout]{*setPrintorder}()}.}
  \item{geneSpotMap}{.}
  \item{plate}{.}
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
  \tab \code{ngrid.r} \tab The number of rows of grids per slide. \cr
  \tab \code{ngrid.c} \tab The number of columns of grids per slide. \cr
  \tab \code{nspot.r} \tab The number of rows of spots in each grid. \cr
  \tab \code{nspot.c} \tab The number of columns of spots in each grid. \cr
  \tab \code{name}    \tab A \code{\link[base]{vector}} of strings which specified the name of each gene. \cr
  \tab \code{id}      \tab A \code{\link[base]{vector}} of strings which specified the id of each gene. \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{anonymize} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{as.data.frame} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.Layout} \tab  -\cr
 \tab \code{\link[aroma:equals.Layout]{equals}} \tab Checks if a Layout object is equals to some other object.\cr
 \tab \code{fromDataFrame} \tab  -\cr
 \tab \code{getBlank} \tab  -\cr
 \tab \code{getGeneGroups} \tab  -\cr
 \tab \code{\link[aroma:getId.Layout]{getId}} \tab Gets the id of one or more spots.\cr
 \tab \code{\link[aroma:getIndex.Layout]{getIndex}} \tab Gets the index of a spot given its location.\cr
 \tab \code{\link[aroma:getIndices.Layout]{getIndices}} \tab Gets the indices of the spots at the given locations.\cr
 \tab \code{getLayoutGroupsByName} \tab  -\cr
 \tab \code{\link[aroma:getLocation.Layout]{getLocation}} \tab Gets the location of a spot given its index.\cr
 \tab \code{\link[aroma:getName.Layout]{getName}} \tab Gets the name of one or more spots.\cr
 \tab \code{getPlate} \tab  -\cr
 \tab \code{getPlateGroups} \tab  -\cr
 \tab \code{getPlateNumber} \tab  -\cr
 \tab \code{\link[aroma:getPosition.Layout]{getPosition}} \tab Gets the position of a set of spots given their indices.\cr
 \tab \code{getPrintdipGroups} \tab  -\cr
 \tab \code{\link[aroma:getPrintorder.Layout]{getPrintorder}} \tab Gets the order of the spots in which they were printed.\cr
 \tab \code{getPrinttipGroups} \tab  -\cr
 \tab \code{getSlideColumnGroups} \tab  -\cr
 \tab \code{getSlideRowGroups} \tab  -\cr
 \tab \code{\link[aroma:gridSize.Layout]{gridSize}} \tab Gets the number of spots in each grid.\cr
 \tab \code{hasIds} \tab  -\cr
 \tab \code{hasNames} \tab  -\cr
 \tab \code{hasPlates} \tab  -\cr
 \tab \code{\link[aroma:indexOf.Layout]{indexOf}} \tab Gets the index of one or more spots from their name or id.\cr
 \tab \code{\link[aroma:nbrOfColumns.Layout]{nbrOfColumns}} \tab Gets the number of columns on a microarray.\cr
 \tab \code{\link[aroma:nbrOfGrids.Layout]{nbrOfGrids}} \tab Gets the number of grids on a microarray.\cr
 \tab \code{nbrOfPlates} \tab  -\cr
 \tab \code{nbrOfReplicates} \tab  -\cr
 \tab \code{\link[aroma:nbrOfRows.Layout]{nbrOfRows}} \tab Gets the number of rows on a microarray.\cr
 \tab \code{\link[aroma:nbrOfSpots.Layout]{nbrOfSpots}} \tab Gets the size of a microarray.\cr
 \tab \code{put} \tab  -\cr
 \tab \code{\link[aroma:read.Layout]{read}} \tab Reads layout information from a tab-delimited file.\cr
 \tab \code{rename} \tab  -\cr
 \tab \code{\link[aroma:set.Layout]{set}} \tab Sets the layout.\cr
 \tab \code{setGeneGroups} \tab  -\cr
 \tab \code{\link[aroma:setId.Layout]{setId}} \tab Sets the id of one or more spots.\cr
 \tab \code{\link[aroma:setName.Layout]{setName}} \tab Sets the name of one or more spots.\cr
 \tab \code{setPlate} \tab  -\cr
 \tab \code{setPlateGroups} \tab  -\cr
 \tab \code{setPrintdipGroups} \tab  -\cr
 \tab \code{setPrintorder} \tab  -\cr
 \tab \code{setPrinttipGroups} \tab  -\cr
 \tab \code{setSlideColumnGroups} \tab  -\cr
 \tab \code{setSlideRowGroups} \tab  -\cr
 \tab \code{\link[aroma:size.Layout]{size}} \tab Gets the size of a microarray.\cr
 \tab \code{\link[aroma:toXYMatrix.Layout]{toXYMatrix}} \tab Layouts out values in the same way as the spots are layout.\cr
 \tab \code{\link[aroma:write.Layout]{write}} \tab Writes the layout information to a tab-delimited file.\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\note{
 There are several functions that returns a Layout object.
 It is only in very special cases that you have to create one yourself.

 In the sma package some functions are related to this class. This,
 class might be backward compatible with these functions, but the reverse
 is not true. The following functions are known be related to this class:
 \code{\link[sma]{init.grid}}, \code{\link[sma]{id2image}} (and \code{image2id}).
}

\details{
  GenBank Accession numbers, SwissProt/TrEMBL Accession numbers or
  Entry Names.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)

  SMA$loadData("mouse.setup")
  layout <- Layout(mouse.setup)
  # or, equivalent...
  layout <- as.Layout(mouse.setup)
}

\seealso{
 \code{\link{MicroarrayData.getLayout}}.
}
\keyword{classes}
