#########################################################################/**
# @RdocClass Replicates
#
# @title "The Replicates class"
#
# \description{
#  @classhierarchy
#
#  The Replicates class provides several useful methods for working 
#  replicated spots on an array.
# }
#
# @synopsis
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#*/#########################################################################
setConstructorS3("Replicates", function(replicates=NULL, type=NA) {
  if (is.matrix(replicates))
    replicates <- matrixToList(replicates)
  else if (!is.list(replicates))
    replicates <- as.list(replicates)

  extend(Object(), "Replicates", 
    replicates = replicates,
    type = type,
    parameter = list()
  )
}, deprecated=TRUE)


setMethodS3("as.character", "Replicates", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste("Replicates: ", sep="");
  s <- paste(s, nbrOfGenes(this), " genes in ", nbrOfSpots(this), " spots", sep="");
  s <- paste(s, " (", formatC(nbrOfSpots(this)/nbrOfGenes(this), width=1, digits=3), " spots/gene)", sep="");
  paste(s, ".", sep="");
}, deprecated=TRUE)




setMethodS3("str", "Replicates", function(object, ...) {
  # To please R CMD check...
  this <- object;

  str(paste(sep="", "a list of ", length(this$replicates), " replicates."), ...);
}, deprecated=TRUE)



setMethodS3("equals", "Replicates", function(this, obj) {
  if (any(nbrOfGenes(this) != nbrOfGenes(obj)))
    return(FALSE);
  if (any(nbrOfSpots(this) != nbrOfSpots(obj)))
    return(FALSE);
  if (any(unlist(this$replicates) != unlist(obj$replicates)))
    return(FALSE);
  TRUE;
}, deprecated=TRUE)


setMethodS3("nbrOfGenes", "Replicates", function(this) {
  length(this$replicates);
}, deprecated=TRUE)


setMethodS3("nbrOfSpots", "Replicates", function(this) {
  length(unlist(this$replicates));
}, deprecated=TRUE)




#########################################################################/**
# @RdocMethod nbrOfReplicates
#
# @title "Gets the number of replicates for each gene"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \value{
#   Returns a @vector of @integers with the number of replicates for each
#   gene.
# }
#
# @author
#
# \seealso{
#  @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfReplicates", "Replicates", function(this, genes=NULL) {
  nbr <- unlist(lapply(this$replicates, FUN=length));
  if (!is.null(genes))
    nbr <- nbr[genes];
  nbr;
}, deprecated=TRUE)




#########################################################################/**
# @RdocMethod hasReplicates
#
# @title "Checks if the microarray(s) has replicates of genes within the slides"
#
# @synopsis
#
# \description{
#   @get "title".
#   Note, for this to work the replicates must be set by using the
#   @seemethod "set" method.
# }
#
# \value{
#   Returns @TRUE if there are replicates of genes, otherwise @FALSE.
# }
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(hasReplicates(layout))     # FALSE
#
#   neighborReplicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
#   set(layout, neighborReplicates)
#
#   print(hasReplicates(layout))     # TRUE
# }
#
# @author
#
# \seealso{
#  To get or set the replicates for one or more genes, see 
#  @seemethod "getReplicates" and @seemethod "set". 
#  To automatically identify replicated spots using the ID's or the Name's,
#  see @seemethod "findReplicates".
#  @seeclass
# }
#*/#########################################################################
setMethodS3("hasReplicates", "Replicates", function(this) {
  (max(nbrOfReplicates(this)) > 1);
}, deprecated=TRUE)





#########################################################################/**
# @RdocMethod getSpot
#
# @title "Gets the spot indices for all replicates of one or more genes"
#
# @synopsis
#
# \arguments{
#   \item{genes}{The genes whose replicates are retrieved. If @NULL, 
#     all genes are considered.}
#   \item{replicates}{The indices within each replicate group of spots
#     to be considered. If @NULL, all spots within each replicate
#     group are considered.}
# }
#
# \value{
#   Returns a @list of length \code{length(genes)} if \code{genes} is
#   given, otherwise all genes are returned. Each component in the @list
#   contains the spot indices for all replicates of that gene. 
# }
#
# \description{
#   @get "title".
# }
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#
#   # Create a matrix where each row represent a gene and the columns
#   # represents two neighboring spot indices. Even if set
#   # expect a list similar to the one returned by getReplicate(), it
#   # also accepts matrices which internally are transformed to lists.
#   neighborReplicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
#   layout$set(neighborReplicates)
#
#   print(getReplicates(layout, 1:4))
# }
#
# @author
#
# \seealso{
#  To check if the slides have replicates or not see 
#  @seemethod "hasReplicates". To set the replicates for one or more genes,
#  see @seemethod "set". To automatically identify replicated
#  spots using the ID's or the Name's, see @seemethod "findReplicates".
#  @seeclass
# }
#*/#########################################################################
setMethodS3("getSpot", "Replicates", function(this, genes=NULL, replicates=NULL) {
  spots <- this$replicates;
  if (is.null(spots))
    throw("No replicates has been specified.");

  if (!is.null(genes))
    spots <- spots[genes];
  
  if (!is.null(spots)) {
    spots <- unlist(lapply(spots, 
      FUN=function(x) {
        vec <- rep(NA, length(replicates));
        replicates <- union(replicates, 1:length(x));
        vec[replicates] <- x[replicates];
        vec;
      }
    ));
  }

  spots;
}, deprecated=TRUE)

setMethodS3("getGene", "Replicates", function(this, spots=NULL) {
  replicates <- this$replicates;
  if (is.null(replicates))
    throw("No replicates has been specified.");

  if (is.null(spots))
    spots <- unlist(replicates)

  reps <- t(listToMatrix(replicates));
  genes <- match(spots, reps);
  (genes-1) %/% nrow(reps) + 1;
}, deprecated=TRUE)




#########################################################################/**
# @RdocMethod set
#
# @title "Sets the replicates for one or more genes"
#
# @synopsis
#
# \description{
#   @get "title". The argument \code{replicates}
#   can either be a list, a matrix or a string with a special value. 
#   If it is a matrix it is first tranformed into a list with all NA's
#   removed. Accepted special strings are \code{"neighboring-pairs"} and
#   \code{"top-bottom"}.
#   For automatically getting the list replicates from the Name's or the
#   ID's the method findReplicates can be used. For more information see
#   @seemethod "findReplicates".
# }
#
# \arguments{
#   \item{replicates}{A @list, a @matrix or a special @character string
#     describing the replicates. The first special string is
#     \code{"neighboring-pairs"} 
#     which specifies that neighboring pairs, (1,2), (3,4), ..., are 
#     replicates. Another special string is \code{"top-bottom"} which
#     specifies that the second half of spots are replicates of the first
#     half of spots.}
#   \item{genes}{The genes whose replicates are set. If @NULL, all
#     genes are considered.}
# }
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#
#   neighborReplicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
#   set(layout, neighborReplicates)
#
#   # Or, identically...
#   set(layout, "neighboring-pairs")
#
#   print(getReplicates(layout, 1:4))
# }
#
# @author
#
# \seealso{
#  To check if the slides have replicates or not see 
#  @seemethod "hasReplicates". To get the replicates for one or more genes, 
#  see @seemethod "getReplicates". To automatically identify replicated
#  spots using the ID's or the Name's, see @seemethod "findReplicates".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("set", "Replicates", function(this, replicates, genes=NULL) {
  # First, make sure 'replicates' is a list and save space by removing NA's.
  if (is.matrix(replicates))
    replicates <- matrixToList(replicates, na.rm=TRUE)
  else if (is.character(replicates)) {
    if (!is.null(genes))
      throw("Argument 'genes' can not be specified if argument 'replicates' is a string.");
    if (replicates == "neighboring-pairs") {
      replicates <- matrix(1:nbrOfSpots(this), ncol=2, byrow=TRUE);
    } else if (replicates == "top-bottom") {
      replicates <- matrix(1:nbrOfSpots(this), ncol=2);
    } else {
      throw("Argument 'replicates' has an unknown value: ", replicates);
    }
    replicates <- matrixToList(replicates, na.rm=TRUE);
  } else if (!is.list(replicates))
    throw("Argument 'replicates' is of an unknown class. Only lists, matrices and certain strings are accepted: ",  data.class(replicates));

  if (!is.null(genes)) {
    reps <- replicates;
    replicates <- this$replicates;
    replicates[genes] <- reps;
  }

  this$replicates <- replicates;
  invisible(this);
}, deprecated=TRUE)



#########################################################################/**
# @RdocMethod fromType
#
# @title "Creates a Replicates object given a type and a layout"
#
# @synopsis
#
# \description{
#   @get "title".
#   The \code{type} argument is a string specifies one of many standard
#   layouts of replicates on a slide.
#   Currently, accepted types are \code{"neighboring-pairs"} and
#   \code{"top-bottom"}.
#   For automatically getting the replicates from the Name's or the
#   ID's see @seemethod "fromLayout".
# }
#
# \arguments{
#   \item{type}{A @character string specifies one of many standard layouts of
#      replicates on a slide. If \code{"neighboring-pairs"}, neighboring
#      pairs, i.e. spots with indices (1,2), (3,4), ..., are set as
#      replicates. If \code{"top-bottom"}, the second half of the spots
#      are replicates of the first half of spots.}
#   \item{layout}{A @see "Layout" object.}
# }
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#
#   # Define neighboring spots to be replicates of the same gene
#   reps1 <- Replicates$fromType("neighboring-pairs", layout)
#
#   # Or, equivalent...
#   neighborReplicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
#   reps2 <- Replicates(neighborReplicates)
#
#   print(equals(reps1, reps2))  # TRUE
# }
#
# @author
#
# \seealso{
#   To automatically identify replicates from the ID's or the Name's of
#   spots (specified by a Layout object), see @seemethod "fromLayout".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("fromType", "Replicates", function(this, type, layout) {
  if (!inherits(layout, "Layout"))
    throw("The argument 'layout' must be of type Layout: ", data.class(layout));

  knownTypes <- c("neighboring-pairs", "top-bottom");
  if (!is.element(type, knownTypes)) {
    throw("The value of 'type' is unknown. Known types are ", knownTypes, ": ", type);
  }

  if (type == "neighboring-pairs") {
    replicates <- matrix(1:nbrOfSpots(layout), ncol=2, byrow=TRUE);
  } else if (type == "top-bottom") {
    replicates <- matrix(1:nbrOfSpots(layout), ncol=2);
  }

  replicates <- matrixToList(replicates, na.rm=TRUE);
  Replicates(replicates, type=type)
}, static=TRUE, deprecated=TRUE);


#########################################################################/**
# @RdocMethod fromLayout
#
# @title "Creates a Replicates object from unique ID's or Name's in a Layout object"
#
# @synopsis
#
# \arguments{
#   \item{layout}{The layout containing the ID's or the Name's.}
#   \item{field}{If \code{"ID"} the ID's are used for identifying the 
#     replicated spots and if \code{"Name"} the Name's are used.}
# }
#
# \description{
#   @get "title".
#   in the given Layout object. The ordering of the identified unique
#   genes will be ordered by their unique spot indices, i.e. as they 
#   appear on the slide.
# }
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"));
#   layout <- getLayout(gpr);
#
#   replicates <- Replicates$fromLayout(layout, "ID")
# }
#
# @author
#
# \seealso{
#   To set the replicates directly in a Layout object see
#   @see "Layout.set".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("fromLayout", "Replicates", function(this, layout, field="ID") {
  if (field == "ID")
    x <- getID(layout)
  else if (field == "Name")
    x <- getName(layout)
  else
    throw("Unknown value of argument 'field'. Valid values are \"ID\" and \"Name\": ", field);

  x.match <- match(x, x);
  replicates <- list();
  for (k in 1:length(x.match)) {
    gene <- x.match[k];
    if (length(replicates) < gene)
      replicates[[gene]] <- k
    else
      replicates[[gene]] <- c(replicates[[gene]], k);
  }

  # Some "genes" may be empty, so remove these.
  empty <- which(lapply(replicates, FUN=length) == 0);
  replicates <- replicates[-empty];

  Replicates(replicates)
}, static=TRUE, deprecated=TRUE);



setMethodS3("setParameter", "Replicates", function(this, paramName, value) {
  replicates <- this$replicates;
  value <- rep(value, length.out=length(replicates));
  values <- list();
  for (k in 1:length(replicates))
    values[[k]] <- rep(value[k], length(replicates[[k]]));
  parameter <- this$parameter;
  if (is.null(parameter)) parameter <- list();
  parameter[[paramName]] <- values;
  this$parameter <- parameter;
  invisible(this);
}, deprecated=TRUE);


setMethodS3("getParameter", "Replicates", function(this, paramName, genes=NULL) {
  value <- this$parameter[[paramName]];
  if (!is.null(genes))
    value <- value[genes];
  value;
}, deprecated=TRUE);



setMethodS3("highlight", "Replicates", function(this, genes=NULL, cex="self", col="self", pch="self", ...) {
  lastPlot <- Device$getPlotParameters();
  object <- lastPlot$object;

  incl <- getSpot(this, genes);

  param <- list();
  for (paramName in c("cex", "col", "pch")) {
    arg <- get(paramName);
    if (length(arg) > 1) {
      arg <- arg[incl];
    } else if (!is.null(arg) && arg == "self") {
      arg <- getParameter(this, paramName, genes);
      arg <- unlist(arg);
    }

    if (length(arg) == 1 && is.na(arg)) arg <- NULL;
#      cat(paramName, "=", arg, "\n");
    param[[paramName]] <- arg;
  } # for (paramName...)

  cex <- param$cex;
  col <- param$col;
  pch <- param$pch;

  highlight(object, include=incl, cex=cex, col=col, pch=pch, ...);
}, deprecated=TRUE);


setMethodS3("text", "Replicates", function(x, genes=NULL, labels=NULL, cex="self", col="self", ...) {
  # To please R CMD check...
  this <- x;

  lastPlot <- Device$getPlotParameters();
  object <- lastPlot$object;

  incl <- getSpot(this, genes);

  param <- list();
  for (paramName in c("cex", "col")) {
    arg <- get(paramName);
    if (length(arg) > 1) {
      arg <- arg[incl];
    } else if (!is.null(arg) && arg == "self") {
      arg <- getParameter(this, paramName, genes);
      arg <- unlist(arg);
    }

    if (length(arg) == 1 && is.na(arg)) arg <- NULL;
#      cat(paramName, "=", arg, "\n");
    param[[paramName]] <- arg;
  } # for (paramName...)

  cex <- param$cex;
  col <- param$col;

  if (!is.null(labels)) {
    if (length(labels) > 1) {
      n <- nbrOfReplicates(this, genes);
      labels <- rep(labels, n);
    }
  }

  text(object, include=incl, labels=labels, cex=cex, col=col, ...);
}, deprecated=TRUE)



############################################################################
# HISTORY:
# 2003-05-04
# o Made all methods deprecated.
# o BUG FIX: The fromLayout() example was trying to read unexisting data.
# 2002-06-24
# * Made the class deprecated.
# 2002-04-21
# * Replaced a few throw()'s with throw()'s.
# 2002-02-27
# * Updated code to make use of setMethodS3's.
# 2001-11-18
# * Added the field 'type'.
# 2001-07-18
# * Created!
############################################################################
