###########################################################################/**
# @RdocClass UnitNamesFile
#
# @title "The UnitNamesFile class"
#
# \description{
#  @classhierarchy
#
#  A UnitNamesFile provides methods for querying the unit names of
#  a given chip type.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "Interface".}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("UnitNamesFile", function(...) {
  extend(Interface(), "UnitNamesFile");
})


setMethodS3("getUnitNames", "UnitNamesFile", abstract=TRUE);

setMethodS3("nbrOfUnits", "UnitNamesFile", function(this, ...) {
  length(getUnitNames(this));
})


setMethodS3("getChipType", "UnitNamesFile", abstract=TRUE);

setMethodS3("getPlatform", "UnitNamesFile", abstract=TRUE);



###########################################################################/**
# @RdocMethod indexOf
#
# @title "Gets the indices of units by their names"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pattern}{A pattern to be used for identifying unit names of 
#      interest.  If @NULL, no regular expression matching is done.}
#   \item{names}{Names to be match exactly to the unit names.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @vector of @integers in [1,N] where N is the number of units
#  in the underlying annotation chip type file.
# }
#
# @author
#
# \seealso{
#   @seemethod "getUnitNames".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("indexOf", "UnitNamesFile", function(this, pattern=NULL, names=NULL, ...) {
  if (!is.null(names)) {
    idxs <- match(names, getUnitNames(this));
  } else if (!is.null(pattern)) {
    idxs <- grep(pattern, getUnitNames(this));
  } else {
    throw("Either argument 'names' or 'pattern' must be specified.");
  }

  idxs;
})


############################################################################
# HISTORY:
# 2008-07-21
# o Renamed UnitNamesInterface to UnitNamesFile.
# 2008-05-18
# o Created.
############################################################################
