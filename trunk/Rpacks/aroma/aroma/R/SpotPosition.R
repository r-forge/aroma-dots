#########################################################################/**
# @RdocClass SpotPosition
#
# @title "Class representing spot positions"
#
# \description{
#  @classhierarchy
#
#   A SpotPosition object contains physical (spatial) coordinates (in units
#   of pixels) of the spots on a microarray slide. 
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{An NxM @matrix with x positions, where N is the number of
#     spots on each slide and M is the number of slides. Alternatively,
#     another SpotPosition object can be given and then argument \code{y}
#     is ignored.}
#   \item{y}{An NxM @matrix with y positions of the same dimension as
#     argument \code{x}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/#########################################################################
setConstructorS3("SpotPosition", function(x=NULL, y=NULL) {
  if (!missing(x)) {
    if (inherits(x, "SpotPosition")) {
      x <- x$x;
      y <- x$y;
    }

    if (is.null(x) || is.null(y))
      throw("Neither 'x' nor 'y' can be NULL.");

    # Assert that x and y are matrices of the same dimensions
    if (!identical(dim(x), dim(y)))
      throw("Argument 'x' and 'y' must be of the same dimensionality.");

    x <- as.matrix(x);
    y <- as.matrix(y);
  } else {
    x <- y <- NA;
  }

  extend(Object(), "SpotPosition", 
    x = x,
    y = y
  )
})

#########################################################################/**
# @RdocMethod as.data.frame
#
# @title "Gets a data frame represenation of the spot positions"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \value{
#  Returns a @data.frame.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("as.data.frame", "SpotPosition", function(x) {
  # To please R CMD check...
  this <- x;

  data.frame(slide=rep(1:nbrOfSlides(this), each=nbrOfSpots(this)),
             spot=rep(1:nbrOfSpots(this), times=nbrOfSlides(this)),
             x=as.vector(this$x), y=as.vector(this$y));
})


setMethodS3("extract", "SpotPosition", function(this, slides=NULL) {
  df <- as.data.frame(this);
  if (!is.null(slides)) {
    incl <- (df$slide %in% slides);
    df <- df[incl,];
  }
  df;
})


setMethodS3("equals", "SpotPosition", function(this, other) {
  if (!inherits(other, "SpotPosition"))
    return(FALSE);
  
  fields <- sort(getFields(this));
  if (any(fields != sort(getFields(other))))
    return(FALSE);

  for (field in fields) {
    if (any(this[[field]] != other[[field]]))
      return(FALSE);
  }
  
  return(TRUE);
})


#########################################################################/**
# @RdocMethod as.character
#
# @title "Gets a character represenation of the spot positions"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("as.character", "SpotPosition", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste(data.class(this), ":", sep="");
  s <- paste(s, " number of spots: ", nbrOfSpots(this), sep="");
  s <- paste(s, " number of slides: ", nbrOfSlides(this), sep="");
  s;
})




#########################################################################/**
# @RdocMethod nbrOfSlides
#
# @title "Gets the number of slides"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfSpots",
#   @seemethod "size".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfSlides", "SpotPosition", function(this) {
  ncol(this$x)
})


#########################################################################/**
# @RdocMethod nbrOfSpots
#
# @title "Gets the number of spots in each of the slides"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfSlides",
#   @seemethod "size".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfSpots", "SpotPosition", function(this) {
  nrow(this$x)
})


#########################################################################/**
# @RdocMethod size
#
# @title "Gets the number of spots in all slides together"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfSlides",
#   @seemethod "nbrOfSpots".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("size", "SpotPosition", function(this) {
  length(this$x)
})





#########################################################################/**
# @RdocMethod getLeftEdge
#
# @title "Gets the x-coordinate of the left most spot"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getLeftEdge", "SpotPosition", function(this) {
  min(this$x, na.rm=TRUE);
})


#########################################################################/**
# @RdocMethod getRightEdge
#
# @title "Gets the x-coordinate of the right most spot"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getRightEdge", "SpotPosition", function(this) {
  max(this$x, na.rm=TRUE);
})



#########################################################################/**
# @RdocMethod getTopEdge
#
# @title "Gets the y-coordinate of the spot at the very top"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getTopEdge", "SpotPosition", function(this) {
  max(this$y, na.rm=TRUE);
})


#########################################################################/**
# @RdocMethod getBottomEdge
#
# @title "Gets the y-coordinate of the spot at the very bottom"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getBottomEdge", "SpotPosition", function(this) {
  min(this$y, na.rm=TRUE);
})



#########################################################################/**
# @RdocMethod getMaxWidth
#
# @title "Gets the maximum x-distance between two spots"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getMaxWidth", "SpotPosition", function(this) {
  abs(getRightEdge(this) - getLeftEdge(this));
})


#########################################################################/**
# @RdocMethod getMaxHeight
#
# @title "Gets the maximum y-distance between two spots"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getMaxHeight", "SpotPosition", function(this) {
  abs(getTopEdge(this) - getBottomEdge(this));
})

#########################################################################/**
# @RdocMethod getAspectRatio
#
# @title "Gets the aspect ratio of the spots"
#
# @synopsis
#
# \description{
#   @get "title", i.e. the ratio between the maximum height and the maximum
#   width. If all spots on an array are included this is \emph{approximately}
#   the aspect ratio of the array itself. This is useful to make graphs that
#   display true physical positions.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getAspectRatio", "SpotPosition", function(this) {
  getMaxHeight(this) / getMaxWidth(this);
})




#########################################################################/**
# @RdocMethod getDistancesTo
#
# @title "Gets the Euclidean distance from one spot to a set of other spots"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getDistancesTo", "SpotPosition", function(this, index) {
  if (!is.numeric(index) || length(index) != 1)
    throw("Argument 'index' must be a numeric scalar.");
  if (index < 0 || index > nrow(this$x))
    throw("Argument 'index' is out of range: ", index);
  
  x0 <- this$x[index,];
  dx <- t(t(this$x) - x0);
  
  y0 <- this$y[index,];
  dy <- t(t(this$y) - y0);
  
  sqrt(dx^2 + dy^2);
})







#########################################################################/**
# @RdocMethod plot
#
# @title "Plots the spot positions in an xy scatter plot"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{xlab, ylab}{@character string or @expression specifying the 
#       x and the y labels.}
#   \item{pch}{The point character to be plotted at each position.}
#   \item{...}{All other arguments accepted by the \code{plot()} function
#       found in the base package.}
# }
#
# @author
#
# \seealso{
#   @seemethod "points"
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plot", "SpotPosition", function(x, xlab="x", ylab="y", pch=20, ...) {
  # To please R CMD check...
  this <- x;

  plot(this$x, -this$y, xlab=xlab, ylab=ylab, pch=pch, ...);
})


#########################################################################/**
# @RdocMethod points
#
# @title "Adds spot positions to an existing plot"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{pch}{The point character to be plotted at each position.}
#   \item{...}{All other arguments accepted by the \code{points()} function
#       found in the base package.}
# }
#
# @author
#
# \seealso{
#   @seemethod "plot"
#   @seeclass
# }
#*/#########################################################################
setMethodS3("points", "SpotPosition", function(x, index=NULL, pch=20, ...) {
  # To please R CMD check...
  this <- x;

  if (is.null(index)) {
    x <-  this$x;
    y <- -this$y;
  } else {
    x <-  this$x[index];
    y <- -this$y[index];
  }
  points(x, y, pch=pch, ...);
})




#########################################################################/**
# @RdocMethod write
#
# @title "Write a SpotPosition object to file"
#
# \description{
#   @get "title". By default, if not overridden
#   by a method in a subclass, it writes the data to a tab-delimited file.
#   Note that subclasses like GenePixData, ScanAlyzeData and SpotData do
#   write files in their special file formats. To force such object of such
#   classes to be written as tab-delimited file, do
#   \code{write.MicroarrayData(object, ...)} instead.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the GPR file to be written.}
#   \item{path}{The path to the GPR file.}
#   \item{overwrite}{If @TRUE, an existing file is overwritten.
#     Otherwise an exception is thrown.}
#   \item{row.names}{If @TRUE, row names are written, otherwise not.}
#   \item{sep}{The separator between the cells.}
#   \item{...}{Other arguments accepted by subclasses or which are passed
#     to \code{write.table}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \seealso{
#   @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "SpotPosition", function(this, filename, path=NULL, slides=NULL, overwrite=FALSE, row.names=FALSE, sep="\t", ..., verbose=FALSE) {
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  # Extracts the fields to be written
  df <- extract(this, slides=slides);
  
  write.table(df, file=filename, row.names=row.names, sep=sep, ...);
});



setMethodS3("read", "SpotPosition", function(this, filename, path=NULL, verbose=FALSE) {
  fields <- c("slide", "spot", "x", "y");
  res <- MicroarrayData$readToList(filename, path=path,
                                   reqFields=fields, verbose=verbose);
  SpotPosition(x=res$x, y=res$y);
}, static=TRUE, trial=TRUE);

 

############################################################################
# HISTORY:
# 2005-10-21
# o Replace 'overwrite' arguments with 'mustNotExist' in calls to Arguments. 
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2004-05-10
# o BUG FIX: SpotPosition$read() was incorrect. Most things I tried to do
#   was already taken care of by MicroarrayData$readToList().
# 2004-05-02
# o Added write() and read().
# o Added as.data.frame() and extract().
# o Added equals().
# 2004-04-26
# o Added getLeftEdge(), ..., getMaxWidth(), ..., getAspectRatio().
# 2004-04-19
# o Added getDistancesTo().
# o Added argument 'index' to points.
# 2003-05-03
# o Added Rdoc comments.
# 2002-04-05
# o Renamed from SpotPositions to SpotPosition.
# 2002-03-24
# o Created.
############################################################################
