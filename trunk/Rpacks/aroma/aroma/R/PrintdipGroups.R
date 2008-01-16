#########################################################################/**
# @RdocClass PrintdipGroups
#
# @title "The PrintdipGroups class"
#
# \description{
#  @classhierarchy
#
#  Do not mistake this class for the @see "PrinttipGroups".
# }
#
# @synopsis
#
# \arguments{
#  \item{layout}{A @see "Layout" object specifying the layout of a set of
#   microarrays.}
#  \item{groups}{A @list of length equal to the number of printdips and
#   each element contains spot indices corresponding to each printdip. If
#   @NULL, the printdips are estimated/guessed from the Layout object.}
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
setConstructorS3("PrintdipGroups", function(layout=NULL, groups=NULL) {
  if (!is.null(layout) && is.null(groups)) {
    groups <- as.list(1:gridSize(layout));
    names(groups) <- 1:gridSize(layout);
  }

  extend(LayoutGroups(layout), "PrintdipGroups",
    groups=groups
  )
})

setMethodS3("as.character", "PrintdipGroups", function(x, ...) {
  # To please R CMD check
  object <- x;

  s <- NextMethod("as.character");
  s <- paste(s, " The size of each printdip groups is ", getSizes(object)[1], ".", sep="");
  s;
})

setMethodS3("getNames", "PrintdipGroups", function(object) {
  names(object$groups);
})

setMethodS3("nbrOfGroups", "PrintdipGroups", function(object) {
  length(object$groups)
})

setMethodS3("getSizes", "PrintdipGroups", function(object) {
  rep(nbrOfGrids(object$layout), nbrOfGroups(object))
})


setMethodS3("getSpots", "PrintdipGroups", function(this, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) {
    groups <- seq(nbrOfGroups(this))
  } else if (is.numeric(groups)) {
    if (any(groups < 1 || groups > nbrOfGroups(this)))
      throw("Argument 'groups' contain a value that is out of range.")
  } else if (is.character(groups)) {
    match <- match(unique(groups), getNames(this))
    ok <- !is.na(match)
    if (sum(!ok) != 0)
      warning("Some of the names asked for where not found among the printdip names.");
    groups <- match[ok];
  }

  dips <- toPrintorderMatrix(this$layout);
  colnames(dips) <- 1:gridSize(this$layout);
  # rownames(dips) <- 1:nbrOfGrids(this$layout); # unnecessary
  dips <- dips[,groups];
  df <- as.data.frame(dips);
  l <- as.list(df);
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})


############################################################################
# HISTORY:
# 2003-09-19
# o Extracted PrintdipGroups into its own source file.
############################################################################
