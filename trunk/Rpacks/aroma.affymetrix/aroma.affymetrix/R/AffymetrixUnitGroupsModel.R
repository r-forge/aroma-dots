###########################################################################/**
# @RdocClass AffymetrixUnitGroupsModel
#
# @title "The AffymetrixUnitGroupsModel class"
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
#   \item{...}{Not used.}
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
setConstructorS3("AffymetrixUnitGroupsModel", function(dataSet=NULL, name="modelUnitGroups", path=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (is.null(path)) {
    if (!is.null(dataSet)) {
      # Path structure: <data-set name>/<model name>/<chip type>/
      # Compare with  : <data-set name>/chip_data/<chip type>/
      path <- filePath(getName(dataSet), name, getChipType(getCdf(dataSet)));
    }
  } 
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path);
  }

  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixFileSet"))
      throw("Argument 'dataSet' is not an AffymetrixFileSet object: ", class(dataSet));
  }

  # Clone the data set.  Since the CDF is defined in the data set,
  # models using different CDF structures must not share the same
  # CDF, thus not the same data set. /HB 2006-08-28
  if (!is.null(dataSet)) {
    dataSet <- clone(dataSet);
  }

  extend(Object(), "AffymetrixUnitGroupsModel",
    path = path,
    dataSet = dataSet,
    parSet = NULL
  )
})

setMethodS3("as.character", "AffymetrixUnitGroupsModel", function(this, ...) {
  s <- NextMethod("as.character", this);
  s <- c(s, sprintf("Dataset: %s.", paste(as.character(getDataSet(this)), collapse=". ")));
  s <- c(s, sprintf("Parameters: %s.", paste(as.character(getParameterSet(this)), collapse=". ")));
  s <- paste(s, collapse=" ");
  s;
})

setMethodS3("getName", "AffymetrixUnitGroupsModel", function(this, ...) {
  # Name from pathname structure: <data set>/<name>/<chip type>/
  path <- getPath(this);
  
  # <data set>/<name>/
  path <- dirname(path);

  # <name>
  name <- basename(path);

  name;
})

setMethodS3("getLabel", "AffymetrixUnitGroupsModel", function(this, ...) {
  label <- this$.label;
  if (is.null(label))
    label <- getName(this, ...);
  label;
})

setMethodS3("setLabel", "AffymetrixUnitGroupsModel", function(this, label, ...) {
  oldLabel <- this$.label;
  this$.label <- label;
  invisible(oldLabel);
})

###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) for this model"
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
#   Returns a @character.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "AffymetrixUnitGroupsModel", function(this, ...) {
  this$path;
})


###########################################################################/**
# @RdocMethod getDataSet
#
# @title "Gets the data set for this model"
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
#   Returns an @see "AffymetrixFileSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDataSet", "AffymetrixUnitGroupsModel", function(this, ...) {
  this$dataSet;
})


###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "AffymetrixUnitGroupsModel", function(this, ...) {
  getCdf(this$dataSet);
})



###########################################################################/**
# @RdocMethod fit
#
# @title "Estimates the model parameters"
#
# \description{
#  @get "title" for the given data set to all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments used by the method for a subclass.}
# }
#
# \value{
#  Returns ???.
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
setMethodS3("fit", "AffymetrixUnitGroupsModel", abstract=TRUE);



setMethodS3("getParameterSet", "AffymetrixUnitGroupsModel", function(this, ...) {
  this$parSet;
})

setMethodS3("setup", "AffymetrixUnitGroupsModel", abstract=TRUE);


############################################################################
# HISTORY:
# 2006-08-28
# o Added getLabel(), which defaults to getName(), and setLabel().
# 2006-08-24
# o Added some Rdoc comments.
# 2006-08-17
# o Created.
############################################################################
