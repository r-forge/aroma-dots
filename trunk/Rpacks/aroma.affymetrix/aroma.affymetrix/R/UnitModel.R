###########################################################################/**
# @RdocClass UnitModel
#
# @title "The UnitModel class"
#
# \description{
#  @classhierarchy
#
#  This class is abstract and represents a generic unit model, i.e.
#  a model that is applied to each unit separately.
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSet}{An @see "AffymetrixCelSet" to which this model should
#      be fitted.}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
#   \item{...}{Arguments passed to the constructor of @see "Model".}
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
setConstructorS3("UnitModel", function(dataSet=NULL, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet"))
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ",
                                                           class(dataSet));
  }

  this <- extend(Model(dataSet=dataSet, ...), "UnitModel");

  # Interpret and append tags
  setTags(this, tags);

  this;
}, abstract=TRUE)



setMethodS3("getTags", "UnitModel", function(this, ...) {
  tags <- this$.tags;
  tags[tags == "*"] <- paste(getAsteriskTag(this), collapse=",");
  tags;
})


setMethodS3("getAsteriskTag", "UnitModel", function(this, ...) {
  "UM";
})


###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Gets the cell indices unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probesets). 
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to \code{getCellIndices()} of the 
#     @see "AffymetrixCdfFile" of the input data set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @list structure consisting of CDF cell indices.
# }
#
# \details{
#   By default, this is just a wrapper function calling
#   \code{getCellIndices()} of the @see "AffymetrixCdfFile"
#   of the input data set.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCellIndices", "UnitModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF cell indices
  verbose && enter(verbose, "Identifying CDF cell indices");
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  cells <- getCellIndices(cdf, ...);
  verbose && exit(verbose);
  
  cells;
})

###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probesets). 
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Arguments passed to \code{getCellIndices()} of the 
#     @see "AffymetrixCdfFile" class (if \code{cdf} was not specified),
#     as well as \code{getUnitIntensities()} of the 
#     @see "AffymetrixCelSet" of the input data set.}
#   \item{force}{If @TRUE, cached cell indices as well as data is re-read.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @list structure that \code{readUnits()} of 
#  @see "AffymetrixCelSet" returns.
# }
#
# \details{
#   The output structure is shaped by the @list structure returned
#   by @seemethod "getCellIndices".  The default is to return whatever
#   the CDF returns, but by overriding the latter, what cells are read
#   and what group-structure each unit has can be fully customized
#   by any subclass.
# }
#
# @author
#
# \seealso{
#   @seemethod "getCellIndices".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "UnitModel", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the CDF cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying CDF cell indices");
  cdfUnits <- getCellIndices(this, units=units, ..., force=force);
  verbose && print(verbose, cdfUnits[1]);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the CEL intensities by units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");
  res <- getUnitIntensities(ds, units=cdfUnits, ..., force=force);
  verbose && exit(verbose);
  verbose && str(verbose, res[1]);

  res;
})


###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies non-fitted units"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments specific to any subclass.}
# }
#
# \value{
#  Returns an @integer @vector of unit indices.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("findUnitsTodo", "UnitModel", abstract=TRUE)

setMethodS3("getFitUnitFunction", "UnitModel", abstract=TRUE, private=TRUE)


############################################################################
# HISTORY:
# 2007-01-06
# o Removed never-used setup() method.
# 2007-01-01
# o Created from former UnitGroupsModel with history as follows:
# 2006-11-19
# o Started to modify methods of this class to work similar to the
#   QuantileNormalizer and AllelicCrosstalkCalibrator classes.
# 2006-09-14
# o Not cloning the data set anymore.  Each model is responsible for 
#   tranforming the data structure their way.  The advantage with this
#   approach is that we can cache read data in the data set object.
# 2006-08-28
# o Added getLabel(), which defaults to getName(), and setLabel().
# 2006-08-24
# o Added some Rdoc comments.
# 2006-08-17
# o Created.

############################################################################
