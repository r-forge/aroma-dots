###########################################################################/**
# @set "class=AffymetrixDataSet"
# @RdocMethod transformOn
# @aliasmethod transformOff
#
# @title "Turns transformations on and off for all data files in the data set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("transformOn", "AffymetrixDataSet", function(this, ...) {
  invisible(lapply(this, FUN=transformOn));
}, protected=TRUE)

setMethodS3("transformOff", "AffymetrixDataSet", function(this, ...) {
  invisible(lapply(this, FUN=transformOff));
}, protected=TRUE)


###########################################################################/**
# @RdocMethod doTransform
#
# @title "Sets if transformation of probe signals should be done or not"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{status}{A @logical.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("doTransform", "AffymetrixDataSet", function(this, status=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  status <- Arguments$getLogical(status);


  oldStatus <- as.logical(this$.doTransform);
  this$.doTransform <- status;

  invisible(oldStatus);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod isTransforming
#
# @title "Checks if transformation is turned on for any of the data files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE or @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isTransforming", "AffymetrixDataSet", function(this, ...) {
  any(unlist(lapply(this, FUN=isTransforming)));
}, protected=TRUE)


###########################################################################/**
# @RdocMethod hasTransforms
#
# @title "Checks if any of the data files has transform functions specified"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE or @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("hasTransforms", "AffymetrixDataSet", function(this, ...) {
  any(unlist(lapply(this, FUN=hasTransforms)));
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-05-15
# o Extracted from AffymetrixDataSet.R.
# o Letting readCelUnits() transform signals improves speed substantially.
# o Making use of new multi-array readCelUnits().
# 2006-02-20
# o Created.
############################################################################
