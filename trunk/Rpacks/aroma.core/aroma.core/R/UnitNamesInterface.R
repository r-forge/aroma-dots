###########################################################################/**
# @RdocClass UnitNamesInterface
#
# @title "The UnitNamesInterface class"
#
# \description{
#  @classhierarchy
#
#  A UnitNamesInterface provides methods for querying the unit names of
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
setConstructorS3("UnitNamesInterface", function(...) {
  extend(Interface(), "UnitNamesInterface");
})


setMethodS3("getUnitNames", "UnitNamesInterface", abstract=TRUE);

setMethodS3("nbrOfUnits", "UnitNamesInterface", function(this, ...) {
  length(getUnitNames(this));
})


setMethodS3("getChipType", "UnitNamesInterface", abstract=TRUE);

setMethodS3("getPlatform", "UnitNamesInterface", abstract=TRUE);



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
setMethodS3("indexOf", "UnitNamesInterface", function(this, pattern=NULL, names=NULL, ...) {
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
# 2008-05-18
# o Created.
############################################################################
