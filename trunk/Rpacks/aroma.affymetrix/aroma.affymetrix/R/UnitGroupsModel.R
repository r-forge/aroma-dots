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
#   \item{dataSet}{The dataset to which this model should be fitted.}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
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
setConstructorS3("UnitGroupsModel", function(dataSet=NULL, tags="", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet"))
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ", class(dataSet));
  }

  this <- extend(Object(), "UnitGroupsModel",
    dataSet = dataSet,
    .tags = NULL,
    parSet = NULL
  );

  setTags(this, tags);

  this;
}, abstract=TRUE)



setMethodS3("as.character", "UnitGroupsModel", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  ds <- getDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  tags <- paste(getTags(ds), collapse=",");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Parameters: %s.", paste(as.character(getParameterSet(this)), collapse=". ")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


###########################################################################/**
# @RdocMethod getRootPath
#
# @title "Gets the root path of this model"
#
# \description{
#  @get "title", which is the class name of the model preceeded by the
#  string \code{"model"}.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getRootPath", "UnitGroupsModel", function(this, ...) {
  # Default root path
  paste("model", class(this)[1], sep="");
})



###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output summarized data set"
#
# \description{
#  @get "title", which is the same as the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "UnitGroupsModel", function(this, ...) {
  ds <- getDataSet(this);
  getName(ds);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this model.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "UnitGroupsModel", function(this, ...) {
  ds <- getDataSet(this);
  c(getTags(ds), this$.tags);
})


setMethodS3("setTags", "UnitGroupsModel", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }
  
  this$.tags <- tags;
})


###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output data set"
#
# \description{
#  @get "title", which is the name with comma separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "UnitGroupsModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of this model"
#
# \description{
#  @get "title" where the parameter files are stored.
#  If non-existing, then the directory is created.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "UnitGroupsModel", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
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
setMethodS3("getDataSet", "UnitGroupsModel", function(this, ...) {
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
setMethodS3("getCdf", "UnitGroupsModel", function(this, ...) {
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
setMethodS3("fit", "UnitGroupsModel", abstract=TRUE);



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

setMethodS3("getParameterSet", "UnitGroupsModel", function(this, ...) {
  this$parSet;
})

setMethodS3("setup", "UnitGroupsModel", abstract=TRUE);


############################################################################
# HISTORY:
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
