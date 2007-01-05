###########################################################################/**
# @RdocClass Preprocessing
#
# @title "The Preprocessing class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a preprocessor that fits a transformation
#  function and transform the data set accordingly.
#  A preprocessor has an input data set, which is transformed into an
#  output data set.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{The input data set as an @see "AffymetrixCelSet".}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \details{
#   Subclasses must implement the \code{process()} method.
# }
#
# @author
#*/###########################################################################
setConstructorS3("Preprocessing", function(dataSet=NULL, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet"))
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ", class(dataSet));
  }


  this <- extend(Object(), "Preprocessing", 
    .tags = tags,
    inputDataSet = dataSet,
    "cached:outputDataSet" = NULL
  );

  setTags(this, tags);

  this;
}, abstract=TRUE)



setMethodS3("getAsteriskTag", "Preprocessing", function(this, ...) {
  # Default '*' tag is the abbreviation from upper-case letters only,
  # e.g. "QuantileNormalization" gives "QN".
  tag <- class(this)[1];
  tag <- strsplit(tag, split="")[[1]];
  tagUC <- toupper(tag);
  keep <- (tag == tagUC);
  tag <- tag[keep];
  tag <- paste(tag, collapse="");
  tag <- sprintf("%s", tag);
  tag;
}, protected=TRUE)


setMethodS3("getRootPath", "Preprocessing", function(this, ...) {
  sprintf("pp%s", capitalize(class(this)[1]));
})



setMethodS3("as.character", "Preprocessing", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  ds <- getInputDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  tags <- paste(getTags(ds), collapse=",");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Number of arrays: %d (%.2fMb)", 
                           nbrOfArrays(ds), getFileSize(ds)/1024^2));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  params <- paste(getParametersAsString(this), collapse=", ");
  s <- c(s, sprintf("Algorithm parameters: (%s)", params));
  s <- c(s, sprintf("Output path: %s", getPath(this)));
  s <- c(s, sprintf("Is done: %s", isDone(this)));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output data set"
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
setMethodS3("getName", "Preprocessing", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this preprocessor.
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
setMethodS3("getTags", "Preprocessing", function(this, ...) {
  tags <- this$.tags;

  ds <- getInputDataSet(this);
  tags <- c(getTags(ds), tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTag(this);

  tags;
})


setMethodS3("setTags", "Preprocessing", function(this, tags="*", ...) {
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
setMethodS3("getFullName", "Preprocessing", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})




setMethodS3("getParametersAsString", "Preprocessing", function(this, ...) {
  params <- getParameters(this);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, protected=TRUE)



setMethodS3("getParameters", "Preprocessing", function(this, ...) {
  NULL;
})




###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the output data set"
#
# \description{
#  @get "title".
#  If non-existing, then the directory is created.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
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
setMethodS3("getPath", "Preprocessing", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell$", "", chipType);  # AD HOC? /HB 2006-12-08

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
# @RdocMethod getInputDataSet
#
# @title "Gets the source data set"
#
# \description{
#  @get "title" that is to be (or has been) preprocessed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getInputDataSet", "Preprocessing", function(this, ...) {
  this$inputDataSet;
})



###########################################################################/**
# @RdocMethod getOutputDataSet
#
# @title "Gets the preprocessed data set"
#
# \description{
#  @get "title", if processed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet" or @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getOutputDataSet", "Preprocessing", function(this, ..., force=FALSE) { 
 outputDataSet <- this$outputDataSet;
  if (force || is.null(outputDataSet)) {
    if (isDone(this)) {
      ds <- getInputDataSet(this);
      clazz <- Class$forName(class(ds)[1]);
      outputDataSet <- clazz$fromFiles(path=getPath(this));
      this$outputDataSet <- outputDataSet;
    }
  }
  outputDataSet;
})


setMethodS3("getOutputFiles", "Preprocessing", function(this, ...) {
  outPath <- getPath(this);
  findFiles(pattern="[.](c|C)(e|E)(l|L)$", paths=outPath, firstOnly=FALSE);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod isDone
#
# @title "Checks if the data set is processed or not"
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
#  Returns @TRUE if the data set is processed, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDone", "Preprocessing", function(this, ...) {
  pathnames <- getOutputFiles(this);
  if (length(pathnames) == 0)
    return(FALSE);

  ds <- getInputDataSet(this);  
  if (length(pathnames) != nbrOfArrays(ds))
    return(FALSE);
  
  return(TRUE);
})


###########################################################################/**
# @RdocMethod process
#
# @title "Processes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already processed is re-processed, 
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "Preprocessing", abstract=TRUE);



############################################################################
# HISTORY:
# 2006-12-20
# o Now isDone() returns FALSE if not all output files are there.  Before
#   an exception was thrown.  This modification allows you to for instance
#   remove a quantile normalized output file, and when reprocessing the
#   data set, only that file will be processed. I made this change after
#   one file was corrupted in a large data set and I did not want to have
#   to reprocess the whole data set.
# 2006-12-08
# o Renamed from PreProcessor.
# 2006-12-07
# o Created from QuantileNormalizer.R.
############################################################################
