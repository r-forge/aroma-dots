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
setConstructorS3("UnitModel", function(dataSet=NULL, tags="", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet"))
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ", class(dataSet));
  }

  this <- extend(Object(), "UnitModel",
    .dataSet = dataSet,
    .tags = NULL,
    parSet = NULL
  );

  setTags(this, tags);

  this;
}, abstract=TRUE)



setMethodS3("as.character", "UnitModel", function(this, ...) {
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
setMethodS3("getRootPath", "UnitModel", function(this, ...) {
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
setMethodS3("getName", "UnitModel", function(this, ...) {
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
setMethodS3("getTags", "UnitModel", function(this, ...) {
  ds <- getDataSet(this);
  c(getTags(ds), this$.tags);
})


setMethodS3("setTags", "UnitModel", function(this, tags="*", ...) {
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
setMethodS3("getFullName", "UnitModel", function(this, ...) {
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
setMethodS3("getPath", "UnitModel", function(this, ...) {
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
setMethodS3("getDataSet", "UnitModel", function(this, ...) {
  this$.dataSet;
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
setMethodS3("getCdf", "UnitModel", function(this, ...) {
  getCdf(getDataSet(this));
})

setMethodS3("getCellIndices", "UnitModel", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF cell indices
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  verbose && enter(verbose, "Identifying CDF cell indices");
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
#  @get "title" for all or a subset of units (probeset) 
#  specially structured for this PLM.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Arguments passed to \code{getCellIndices()} of the 
#     @see "AffymetrixCdfFile" class (if \code{cdf} was not specified),
#     but also to the \code{readUnits()} method of the 
#     @see "AffymetrixCelSet" class.}
# }
#
# \value{
#  Returns the @list structure that \code{readUnits()} of 
#  @see "AffymetrixCelSet" returns.
# }
#
# @author
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "UnitModel", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF cell indices
  verbose && enter(verbose, "Identifying CDF cell indices");
  cdfUnits <- getCellIndices(this, units=units, ...);
  verbose && print(verbose, cdfUnits[1]);
  verbose && exit(verbose);

  # Get the CEL intensities by units
  ds <- getDataSet(this);
  verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");
  res <- getUnitIntensities(ds, units=cdfUnits, ...);
  verbose && exit(verbose);
  verbose && str(verbose, res[1]);

  res;
})

setMethodS3("getFitUnitFunction", "UnitModel", abstract=TRUE, protected=TRUE);


setMethodS3("findUnitsTodo", "UnitModel", abstract=TRUE);


setMethodS3("getParameterSet", "UnitModel", function(this, ...) {
  this$parSet;
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
setMethodS3("fit", "UnitModel", abstract=TRUE);



############################################################################
# HISTORY:
# 2007-01-01
# o Created from UnitGroupsModel.R.
############################################################################
