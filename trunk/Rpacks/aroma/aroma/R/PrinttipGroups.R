#########################################################################/**
# @RdocClass PrinttipGroups
#
# @title "The PrinttipGroups class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{layout}{A @see "Layout" object specifying the layout of a set of
#   microarrays.}
#  \item{groups}{A @list of length equal to the number of printtips and
#   each element contains spot indices corresponding to each printtip. If
#   @NULL, the printtips are obtained from the Layout object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# @examples "../incl/LayoutGroups.Rex"
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setConstructorS3("PrinttipGroups", function(layout=NULL, groups=NULL) {
  if (!is.null(layout) && is.null(groups))
    groups <- as.list(1:nbrOfGrids(layout));

  extend(LayoutGroups(layout), "PrinttipGroups",
    groups=groups
  )
})

setMethodS3("as.character", "PrinttipGroups", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  idx <- getSpots(this);
  first <- unlist(lapply(idx, FUN=function(x) x[1]));
  last <- unlist(lapply(idx, FUN=function(x) x[length(x)]));
  tmp <- paste(first, last, sep="-", collapse=", ");
  s <- paste(s, " The printtip groups are ", tmp, ".", sep="");
  s;
})

setMethodS3("getNames", "PrinttipGroups", function(this) {
  names(this$groups);
})

setMethodS3("nbrOfGroups", "PrinttipGroups", function(this) {
  length(this$groups)
})

setMethodS3("getSizes", "PrinttipGroups", function(this) {
  unlist(lapply(this$groups, FUN=length))
})


setMethodS3("getSpots", "PrinttipGroups", function(this, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) {
    groups <- seq(nbrOfGroups(this))
  } else if (is.numeric(groups)) {
    if (any(groups < 1 || groups > nbrOfGroups(this)))
      throw("Argument 'groups' contain a value that is out of range.")
  } else if (is.character(groups)) {
    match <- match(unique(groups), getNames(this))
    ok <- !is.na(match)
    if (sum(!ok) != 0)
      warning("Some of the names asked for where not found among the printtip names.");
    groups <- match[ok];
  }
  
  grids <- unlist(this$groups[groups]);
  gridSize <- gridSize(this$layout);

  gridOffset <- (grids-1) * gridSize;
  idx <- matrix((1:gridSize), nrow=length(grids), ncol=gridSize, byrow=TRUE);
  idx <- idx + gridOffset;
  l <- as.list(as.data.frame(t(idx), optional=TRUE));
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})


############################################################################
# HISTORY:
# 2003-09-19
# o Extracted PrinttipGroups into its own source file.
############################################################################
