###########################################################################/**
# @RdocClass Model
#
# @title "The Model class"
#
# \description{
#  @classhierarchy
#
#  This class is abstract and represents a generic model that applies
#  to a data set.
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSet}{The data set to which this model should be fitted.}
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
setConstructorS3("Model", function(dataSet=NULL, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet"))
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ",
                                                           class(dataSet));
  }

  this <- extend(Object(), "Model",
    .dataSet = dataSet,
    .tags = NULL,
    parSet = NULL
  );

  # Interpret and append tags
  setTags(this, tags);

  this;
}, abstract=TRUE)



setMethodS3("as.character", "Model", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  ds <- getDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  tags <- paste(getTags(ds), collapse=",");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Parameters: %s.", 
                  paste(as.character(getParameterSet(this)), collapse=". ")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


###########################################################################/**
# @RdocMethod getRootPath
#
# @title "Gets the root path of this model"
#
# \description{
#  @get "title".
#  By default, this is the string \code{"model"} appended by the capitalized
#  name of the model class, e.g. \code{"modelUnit"}.
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
setMethodS3("getRootPath", "Model", function(this, ...) {
  # Default root path
  paste("model", class(this)[1], sep="");
})



###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output data set"
#
# \description{
#  @get "title", which is the same as the name of the input data set.
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
setMethodS3("getName", "Model", function(this, ...) {
  ds <- getDataSet(this);
  getName(ds);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which consists of the tags of the input data set followed
#  by an additional set of tags added by the model.
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
#   @seemethod "setTags"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "Model", function(this, ...) {
  ds <- getDataSet(this);
  c(getTags(ds), this$.tags);
})



###########################################################################/**
# @RdocMethod setTags
#
# @title "Sets the tags to be appended"
#
# \description{
#  @get "title" to the tags of the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{tags}{A @character @vector of tags.
#    The tags may also be passed as comma-separated strings.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getTags"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setTags", "Model", function(this, tags=NULL, ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})


###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output set"
#
# \description{
#  @get "title", which consists of the name with appended 
#  comma-separated tags.
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
setMethodS3("getFullName", "Model", function(this, ...) {
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
#  @get "title" where the parameter files are located.
# }
#
# @synopsis
#
# \arguments{
#  \item{mkdirs}{If @TRUE, the directory is created, if missing.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \details{
#  If the path does not exist, it is created.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "Model", function(this, mkdirs=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell$", "", chipType);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Create directory?
  if (mkdirs) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
})


###########################################################################/**
# @RdocMethod getDataSet
#
# @title "Gets the input data set for this model"
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
#   Returns an @see "AffymetrixCelSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDataSet", "Model", function(this, ...) {
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
setMethodS3("getCdf", "Model", function(this, ...) {
  getCdf(getDataSet(this));
}, private=TRUE)


setMethodS3("getParameterSet", "Model", function(this, ...) {
  this$parSet;
}, private=TRUE)


###########################################################################/**
# @RdocMethod fit
#
# @title "Estimates the model parameters"
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
#  Returns an @integer @vector specifying what units where fitted.
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
setMethodS3("fit", "Model", abstract=TRUE);


setMethodS3("getLabel", "Model", function(this, ...) {
  label <- this$.label;
  if (is.null(label))
    label <- getName(this, ...);
  label;
}, private=TRUE)

setMethodS3("setLabel", "Model", function(this, label, ...) {
  oldLabel <- this$.label;
  this$.label <- label;
  invisible(oldLabel);
}, private=TRUE)


############################################################################
# HISTORY:
# 2007-01-06
# o Extracted from UnitModel.R.  Intended to cover all types of models,
#   e.g. RmaPlm, GladModel, PlasqModel etc.
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
