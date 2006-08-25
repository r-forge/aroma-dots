###########################################################################/**
# @set "class=AffymetrixDataSet"
# @RdocMethod getReadMap
#
# @title "Gets the read map for the data set"
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
#   Returns an @integer @vector, if a read map exists, otherwise @NULL.
# }
#
# @author
#
# \seealso{
#   @seemethod "setApdMap" and @seemethod "hasApdMap".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getReadMap", "AffymetrixDataSet", function(this, ...) {
  getReadMap(this$dataFiles[[1]]);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod setApdMap
#
# @title "Sets the APD map for the data set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{apdMap}{An @see "ApdMap" object, or @NULL.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getReadMap" and @seemethod "hasApdMap".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setApdMap", "AffymetrixDataSet", function(this, apdMap, ...) {
  lapply(this$dataFiles, FUN=setApdMap, apdMap);
  invisible(this);
}, protected=TRUE)

setMethodS3("compactApdMap", "AffymetrixDataSet", function(this, ...) {
  apdMap <- getApdMap(this);
  setApdMap(this, apdMap=apdMap, ...);
}, protected=TRUE)

setMethodS3("getApdMap", "AffymetrixDataSet", function(this, ...) {
  getApdMap(this$dataFiles[[1]], ...);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod hasApdMap
#
# @title "Checks if the data set has a read map"
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
#   Returns @TRUE if a read map exists, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "getReadMap" and @seemethod "setApdMap".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("hasApdMap", "AffymetrixDataSet", function(this, ...) {
  hasApdMap(this$dataFiles[[1]]);
}, protected=TRUE)


setMethodS3("createOptimalReadMap", "AffymetrixDataSet", function(this, path...) {
  createOptimalReadMap(this$dataFiles[[1]], ...);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-05-29
# o Added compactApdMap().  Update setApdMap() to set the APD map of all
#   data files, to minimize the risk for filling up the memory.
# 2006-05-15
# o Extracted from AffymetrixDataSet.R.
# 2006-04-09
# o Now the read map is loaded automatically when fromFiles() used.
# 2006-03-04
# o Added mapping functions.
############################################################################
