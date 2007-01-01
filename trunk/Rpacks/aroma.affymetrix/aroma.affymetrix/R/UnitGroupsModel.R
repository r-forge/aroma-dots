###########################################################################/**
# @RdocClass UnitGroupsModel
#
# @title "The UnitGroupsModel class"
#
# \description{
#  @classhierarchy
#
#  This class is abstract and represents a generic unit-group model, i.e.
#  a model that applies to each group in each unit.  For instance,
#  most probeset-summary models such as the RMA model and 
#  the Li & Wong model belongs to this class of models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to constructor of superclass @see "UnitModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
# }
#*/###########################################################################
setConstructorS3("UnitGroupsModel", function(...) {
  extend(UnitModel(...), "UnitGroupsModel")
}, abstract=TRUE)


setMethodS3("getLabel", "UnitGroupsModel", function(this, ...) {
  label <- this$.label;
  if (is.null(label))
    label <- getName(this, ...);
  label;
})

setMethodS3("setLabel", "UnitGroupsModel", function(this, label, ...) {
  oldLabel <- this$.label;
  this$.label <- label;
  invisible(oldLabel);
})

setMethodS3("setup", "UnitGroupsModel", abstract=TRUE);


############################################################################
# HISTORY:
# 2007-01-01
# o Now inheriting from new superclass UnitModel.
# 2006-11-19
# o Started to modify methods of this class to work similar to the
#   QuantileNormalizer and AllelicCrosstalkCalibrator classes.
# 2006-09-14
# o Not cloning the dataset anymore.  Each model is responsible for 
#   tranforming the data structure their way.  The advantage with this
#   approach is that we can cache read data in the dataset object.
# 2006-08-28
# o Added getLabel(), which defaults to getName(), and setLabel().
# 2006-08-24
# o Added some Rdoc comments.
# 2006-08-17
# o Created.
############################################################################
