# SLOW # setMethodS3("getPosition", "Layout", function(this, index=NULL) {
# SLOW #   if (is.null(index)) {
# SLOW #     index <- 1:nbrOfSpots(this);
# SLOW #   } else if (any(index < 1))
# SLOW #     throw("index is out of range.");
# SLOW # 
# SLOW #   n <- nbrOfSpots(this);   # Modulo n just in case
# SLOW #   index <- ((index-1) %% n)+1;
# SLOW # 
# SLOW #   location <- this$getLocation(index);
# SLOW # 
# SLOW #   nspot.r <- this$nspot.r;
# SLOW #   nspot.c <- this$nspot.c;
# SLOW #   positions <- c();
# SLOW #   for (k in 1:nrow(location)) {
# SLOW #     row <- (location[k,1]-1)*nspot.r+location[k,3];
# SLOW #     column <- (location[k,2]-1)*nspot.c+location[k,4];
# SLOW #     position <- c(row, column);
# SLOW #     positions <- c(positions, position);
# SLOW #   }
# SLOW #   matrix(positions, ncol=2, byrow=TRUE);
# SLOW # })
# SLOW # 
# SLOW # 
# SLOW # # This new (2001-07-09) approach is about four times faster!
# SLOW # setMethodS3("getPosition", "Layout", function(this, index=NULL) {
# SLOW #   if (is.null(index)) {
# SLOW #     index <- 1:nbrOfSpots(this);
# SLOW #   } else if (any(index < 1))
# SLOW #     throw("index is out of range.");
# SLOW # 
# SLOW #   n <- nbrOfSpots(this);
# SLOW #   index <- ((index-1) %% n)+1; # Modulo n just in case
# SLOW #   posmap <- toXYMatrix(this, 1:n);
# SLOW #   # The match below is the bottle neck; everything else is FAST!
# SLOW #   idxmap <- match(index, posmap);
# SLOW #   rows <- ((idxmap-1) %% nbrOfRows(layout) + 1);
# SLOW #   cols <- ((idxmap-1) %/% nbrOfRows(layout) + 1);
# SLOW #   matrix(c(rows,cols), ncol=2, byrow=FALSE);
# SLOW # })

# Third approach in Layout.R is about seven times faster than the second
# approach above and 30 times the first approach!


############################################################################
############################################################################
## 
##  DUPLICATES/REPLICATES FUNCTIONS
## 
############################################################################
############################################################################

setMethodS3("nbrOfReplicates", "Layout", function(object) {
  if (!hasReplicates(object)) 
    1
  else {
    rep <- getReplicates(object);
    max(nbrOfReplicates(rep), na.rm=TRUE);
  }
});


#########################################################################/**
# @set "class=Layout"
# @RdocMethod hasReplicates
#
# \title{Checks if the microarray(s) has replicates of genes within the slides}
#
# @synopsis
#
# \description{
#   Checks if the microarray(s) has replicates of genes within the slides.
#   Note, for this to work the replicates must be set by using the
#   setReplicates method (see \link{Layout.setReplicates}).
# }
#
# \value{
#   Returns \@TRUE if there are replicates of genes, otherwise 
#   \@FALSE.
# }
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(hasReplicates(layout))     # FALSE
#
#   neighborReplicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
#   layout$setReplicates(neighborReplicates)
#
#   print(hasReplicates(layout))     # TRUE
# }
#
# @author
#
# \seealso{
#  To get or set the replicates for one or more genes, see 
#  \link{Layout.getReplicates} and \link{Layout.setReplicates}. 
#  To automatically identify replicated spots using the ID's or the Name's,
#  see \link{Layout.findReplicates}.
#  For more methods in this class see \code{\link{Layout}}.
# }
#*/#########################################################################
setMethodS3("hasReplicates", "Layout", function(this) {
  !is.null(this$replicates);
}, deprecated=TRUE)


#########################################################################/**
# @RdocMethod getReplicates
#
# \title{Gets the spot indices for all replicates of one or more genes}
#
# @synopsis
#
# \arguments{
#   \item{genes}{The genes whose replicates are retrieved. If \@NULL, 
#     all genes are considered. Default value is \@NULL.}
# }
#
# \value{
#   Returns a list of length \code{length(genes)} if \code{genes} is
#   given, otherwise all genes are returned. Each component in the list
#   contains the spot indices for all replicates of that gene. 
# }
#
# \description{
#   Gets the spot indices for all replicates of one or more genes.
# }
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#
#   # Create a matrix where each row represent a gene and the columns
#   # represents two neighboring spot indices. Even if setReplicates
#   # expect a list similar to the one returned by getReplicate(), it
#   # also accepts matrices which internally are transformed to lists.
#   neighborReplicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
#   layout$setReplicates(neighborReplicates)
#
#   print(getReplicates(layout, 1:4))
# }
#
# @author
#
# \seealso{
#  To check if the slides have replicates or not see 
#  \link{Layout.hasReplicates}. To set the replicates for one or more genes,
#  see \link{Layout.setReplicates}. To automatically identify replicated
#  spots using the ID's or the Name's, see \link{Layout.findReplicates}.
#  For more methods in this class see \code{\link{Layout}}.
# }
#*/#########################################################################
setMethodS3("getReplicates", "Layout", function(this) {
  this$replicates;
}, deprecated=TRUE)





#########################################################################/**
# @RdocMethod setReplicates
#
# \title{Sets the replicates for one or more genes}
#
# @synopsis
#
# \description{
#   Sets the replicates for one or more genes. The argument \code{replicates}
#   can either be a list, a matrix or a string with a special value. 
#   If it is a matrix it is first tranformed into a list with all NA's
#   removed. Accepted special strings are \code{"neighboring-pairs"} and
#   \code{"top-bottom"}.
#   For automatically geting the list replicates from the Name's or the
#   ID's the method findReplicates can be used. For more information see
#   \link{Layout.findReplicates}.
# }
#
# \arguments{
#   \item{replicates}{A list, a matrix or a special string describing the
#     replicates. The first special string is \code{"neighboring-pairs"} 
#     which specifies that neighboring pairs, (1,2), (3,4), ..., are 
#     replicates. Another special string is \code{"top-bottom"} which
#     specifies that the second half of spots are replicates of the first
#     half of spots.}
#   \item{genes}{The genes whose replicates are set. If \@NULL, all
#     genes are considered. Default value is \@NULL.}
# }
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#
#   neighborReplicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
#   setReplicates(layout, neighborReplicates)
#
#   # Or, identically...
#   setReplicates(layout, "neighboring-pairs")
#
#   print(getReplicates(layout, 1:4))
# }
#
# @author
#
# \seealso{
#  To check if the slides have replicates or not see 
#  \link{Layout.hasReplicates}. To get the replicates for one or more genes, 
#  see \link{Layout.getReplicates}. To automatically identify replicated
#  spots using the ID's or the Name's, see \link{Layout.findReplicates}.
#   For more methods in this class see \code{\link{Layout}}.
# }
#*/#########################################################################
setMethodS3("setReplicates", "Layout", function(this, replicates) {
  if (inherits(replicates, "Replicates")) {
    this$replicates <- replicates;
  } else if (is.character(replicates)) {
    if (replicates == "ID" || replicates == "Name")
      replicates <- fromLayout.Replicates(layout=this, field=replicates)
    else {
      replicates <- fromType.Replicates(type=replicates, layout=this);
      this$geneSpotMap <- replicates;
    }
    this$replicates <- replicates;
  }
  invisible(this);
}, deprecated=TRUE)



setMethodS3("getGeneIndex", "Layout", function(this, spots=NULL) {
  if (!hasReplicates(this))
    throw("Could not get replicated genes. Replicates are not specified for this layout.");

  getGene(getReplicates(this), spots=spots);
}, deprecated=TRUE)



setMethodS3("getSpotIndex", "Layout", function(this, genes=NULL, replicates=NULL) {
  if (!hasReplicates(this))
    throw("Could not get replicated spots. Replicates are not specified for this layout.");

  getSpot(getReplicates(this), genes=genes, replicates=replicates);
}, deprecated=TRUE)




#########################################################################/**
# @RdocMethod getGeneReplicateIndex
#
# \title{Returns the index of the spots as (gene,replicate) indices.}
#
# @synopsis
#
# \description{
#   Returns the index of the spots as (gene,replicate) indices.
# }
#
# \value{
#   Returns a matrix where the rows represents the genes and the columns the
#   replicates.
# }
#
# \details{
# }
#
# \examples{
#   layout <- Layout(2,2, 3,3)
#   setReplicates(layout, "neighboring-pairs")
#   gr <- getGeneReplicateIndex(layout)
#
#   #        [,1] [,2]
#   #   [1,]    1    2
#   #   [2,]    3    4
#   #   [3,]    5    6
#   #   [4,]    7    8
#   #   [5,]    9   10
#   #   [6,]   11   12
#   #   [7,]   13   14
#   #   [8,]   15   16
#   #   [9,]   17   18
#   #  [10,]   19   20
#   #  [11,]   21   22
#   #  [12,]   23   24
#   #  [13,]   25   26
#   #  [14,]   27   28
#   #  [15,]   29   30
#   #  [16,]   31   32
#   #  [17,]   33   34
#   #  [18,]   35   36
# }
#
# @author
#
# \seealso{
#    \code{\link{Microarray.getGeneSlideReplicateIndex}}.
# }
#*/#########################################################################
setMethodS3("getGeneReplicateIndex", "Layout", function(this) {
  rep <- getReplicates(this)$replicates;
  if (is.null(rep))
    throw("Could not get replicated genes because replicates are not set.");
  rep <- listToMatrix(rep);
  rep;
}, private=TRUE, deprecated=TRUE)



############################################################################
# HISTORY:
# 2002-05-03
# o Extracted a lot of methods identified to be obsolete into file
#   Layout.obsolete.R
# 2002-05-02
# o TYPO FIX: Type fix in string returned by as.character().
# 2002-05-01
# * Added field plate, getPlate(), setPlate(), hasPlate(), getPlateNumber()
#   and nbrOfPlates().
# * Added as.data.frame(), read() and write().
# 2002-04-21
# * Added getGeneReplicateIndex() and getGeneSlideReplicateIndex()
#   to MicroarrayData.
# * Removed obsolete getReferenceFields().
# 2002-04-13
# * Added optional arguments ignoreCase=TRUE and regexpr=FALSE to indexOf().
#   Updated the Rdoc accordingly and added a few more examples.
#   This update was trigged by questions from Lei Jiang.
# 2002-04-05
# * Added default values of 'index' in getLocation() and getPosition().
# 2002-03-29
# * BUG FIX: Added as.character() to getID() and getName() just to make sure
#   it is not a factor that is returned. Will fix this in GenePixData.R too.
#   This is a double security so it won't happend again.
# 2002-02-26
# * Updated the Rdoc's.
# 2001-01-24
# * Modified source to make use of setMethodS3 and setClassS3.
# * Added nbrOfReplicates().
# 2001-11-18
# * Freshend the code of setReplicates().
# 2001-08-08
# * Added toPrintorderMatrix().
# 2001-08-06
# * BUG FIX: equals() didn't work correctly.
# 2001-08-01
# * Moving the functionalites of replicates back to Layout. For now only
#   the API, but later also the internal code... maybe.
# 2001-07-18
# * Renamed the fields .name and .id to name and id.
# * Created a new class Replicates and moved all replicate/duplicate 
#   functionalities there.
# 2001-07-15
# * Added findReplicates().
# * Added support to get, set and check replicates of genes. Will have to
#   define the concept of a gene index, but it should be pretty straight-
#   forward. Also, I would like to add a method which automatically
#   generates the replicate matrix by look at the list of Names or IDs.
# 2001-07-12
# * Added get- and setName(). GenePix is using both a Name and a ID field.
# 2001-07-11
# * Updated the Rdoc comments.
# 2001-07-09
# * Improved the speed of getPosition() by 30 times!
# 2001-07-07
# * Moved all duplicates/replicates methods in MicroarrayData here.
# * Made all function take modulo n of all indices.
# 2001-07-04
# * Added toXYMatrix(). Useful for plotSpatial() etc.
# 2001-07-01
# * Updated some of the Rdoc comments.
# 2001-06-30
# * Added id's to the Layout class.
# 2001-06-24
# * Made getLocation() work on vectors too. 
# * Added some Rd comments.
# 2001-05-14
# * Added getInternalReferences() for improving gco() performance.
# 2001-04-08
# * Bug fix: Argument indices in getPositions(indices) was misspelled.
# 2001-04-02
# * Added getIndices(),getLocations(), getPostion(), getPositions().
#   However, these methods need some speed optimization because now they're
#   quite slow. I don't wanna spend more time on this right now!
# 2001-03-25
# * Added Rdoc comments.
# 2001-03-19
# * Created.
############################################################################
