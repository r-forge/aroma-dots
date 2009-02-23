###########################################################################/**
# @RdocClass AromaTransform
#
# @title "The AromaTransform class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a transform (algorithm/operator) that 
#  transforms data.  A transform has an input data set, which is 
#  transformed into an output data set.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{The input data set as an @see "AromaMicroarrayDataSet".}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
#   \item{...}{Not used.}
#   \item{.reqSetClass}{Internal argument.}
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
setConstructorS3("AromaTransform", function(dataSet=NULL, tags="*", ..., .reqSetClass="AromaMicroarrayDataSet") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, .reqSetClass)) {
      throw("Argument 'dataSet' is not an ", .reqSetClass, " object: ", 
                                                          class(dataSet)[1]);
    }
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }


  this <- extend(Object(), "AromaTransform", 
    .tags = tags,
    .inputDataSet = dataSet,
    "cached:.outputDataSet" = NULL
  );

  setTags(this, tags);

  this;
}, abstract=TRUE)


setMethodS3("clearCache", "AromaTransform", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".outputDataSet")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)




setMethodS3("getAsteriskTags", "AromaTransform", function(this,...) {
  # Create a default asterisk tags for any class by extracting all
  # capital letters and pasting them together, e.g. AbcDefGhi => ADG.
  name <- class(this)[1];

  # Remove any 'Model' suffixes
  name <- gsub("Model$", "", name);

  name <- capitalize(name);

  # Vectorize
  name <- strsplit(name, split="")[[1]];

  # Identify upper case
  name <- name[(toupper(name) == name)];

  # Paste
  name <- paste(name, collapse="");

  tags <- name;

  tags;
}, private=TRUE)


setMethodS3("getRootPath", "AromaTransform", function(this, ...) {
  sprintf("pp%s", capitalize(class(this)[1]));
}, private=TRUE)



setMethodS3("as.character", "AromaTransform", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  ds <- getInputDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  tags <- paste(getTags(ds), collapse=",");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("User tags: %s", paste(this$.tags, collapse=",")));
  s <- c(s, sprintf("Asterisk ('*') tags: %s", getAsteriskTags(this, collapse=",")));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Number of files: %d (%.2fMB)", 
                           nbrOfFiles(ds), getFileSize(ds)/1024^2));
  s <- c(s, sprintf("Platform: %s", getPlatform(ds)));
  s <- c(s, sprintf("Chip type: %s", getChipType(ds)));
  params <- paste(getParametersAsString(this), collapse=", ");
  s <- c(s, sprintf("Algorithm parameters: (%s)", params));
  s <- c(s, sprintf("Output path: %s", getPath(this)));
  s <- c(s, sprintf("Is done: %s", isDone(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


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
setMethodS3("getName", "AromaTransform", function(this, ...) {
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
#  of this transformation.
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
setMethodS3("getTags", "AromaTransform", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this);
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "AromaTransform", function(this, tags="*", ...) {
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
setMethodS3("getFullName", "AromaTransform", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



setMethodS3("getParametersAsString", "AromaTransform", function(this, ...) {
  params <- getParameters(this, expand=FALSE);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, private=TRUE)



setMethodS3("getParameters", "AromaTransform", function(this, ...) {
  NULL;
}, private=TRUE)




###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the output data set"
#
# \description{
#  @get "title".
#  If non-existing, then the directory is created.
#  Windows Shortcut links are recognized.
# }
#
# @synopsis
#
# \arguments{
#   \item{create}{If @TRUE, the path is created, otherwise not.}
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
setMethodS3("getPath", "AromaTransform", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this));
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  # Create path?
  if (create) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
})



###########################################################################/**
# @RdocMethod getInputDataSet
#
# @title "Gets the input data set"
#
# \description{
#  @get "title" that is to be (or has been) transformed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AromaMicroarrayDataSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getInputDataSet", "AromaTransform", function(this, ...) {
  this$.inputDataSet;
})




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
#   \item{verbose}{See @see "R.utils::Verbose".}
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
setMethodS3("isDone", "AromaTransform", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Checking if data set is \"done\"");

  pathnames <- getOutputFiles(this);
  if (length(pathnames) == 0) {
    verbose && cat(verbose, "NOT done. No output files found.");
    return(FALSE);
  }

  ds <- getInputDataSet(this);  

  if (length(pathnames) < nbrOfFiles(ds)) {
    verbose && cat(verbose, "NOT done. Too few output files: ", 
                                   length(pathnames), " < ", nbrOfFiles(ds));
    return(FALSE);
  }

  if (length(pathnames) > nbrOfFiles(ds)) {
    throw("Too many output files found: ", 
                                  length(pathnames), " > ", nbrOfFiles(ds));
  }

  verbose && cat(verbose, "Done. All output files are there: ", 
                                                          length(pathnames));

  verbose && exit(verbose);

  return(TRUE);
})



setMethodS3("getOutputFiles", "AromaTransform", function(this, pattern=NULL, ...) {
  # Argument 'pattern':
  if (is.null(pattern)) {
    # Default filename pattern find non-private (no dot prefix) files with
    # the same file name extension as the input data set.
    ds <- getInputDataSet(this);
    df <- getFile(ds, 1);
    pattern <- sprintf("^[^.].*[.]%s$", getFilenameExtension(df));
  } else {
    pattern <- Arguments$getRegularExpression(pattern=pattern);
  }

  outPath <- getPath(this);
  affxparser::findFiles(pattern=pattern, paths=outPath, firstOnly=FALSE);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getOutputDataSet
#
# @title "Gets the transformed data set"
#
# \description{
#  @get "title", if processed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to static method \code{fromFiles()} of
#      the class of the input @see "AromaMicroarrayDataSet".}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "AromaMicroarrayDataSet" or @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("getOutputDataSet", "AromaTransform", function(this, ..., force=FALSE, verbose=FALSE) { 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting output data set for ", class(this)[1]);

  outputDataSet <- this$.outputDataSet;

  if (force || is.null(outputDataSet)) {
    verbose && enter(verbose, "Checking to see if data set is \"done\"");
    isDone <- isDone(this, verbose=less(verbose));
    verbose && exit(verbose);

    if (isDone) {
      verbose && enter(verbose, "Retrieving input data set");
      ds <- getInputDataSet(this);
      verbose && exit(verbose);
      verbose && enter(verbose, "Retrieving files for ", class(ds)[1], 
                                                        " output data set");

      # Inherit the CDF from the input data set.
      args <- list(path=getPath(this), ...);

      clazz <- Class$forName(class(ds)[1]);
      staticMethod <- clazz$fromFiles;
      args$verbose <- less(verbose);
      outputDataSet <- do.call("staticMethod", args=args);
      rm(staticMethod, args); # Not needed anymore

      verbose && exit(verbose);

      this$.outputDataSet <- outputDataSet;
    }
  }
  verbose && exit(verbose);

  outputDataSet;
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
setMethodS3("process", "AromaTransform", abstract=TRUE);





############################################################################
# HISTORY:
# 2008-05-31
# o Updated the default filename pattern for getOutputDataFiles().
# 2008-05-25
# o Added generic getOutputDataSet() and getOutputDataFiles().
# 2008-05-24
# o Added argument 'create' to getPath().
# o Added user tags and asterisk tags to the as.character() output.
# 2008-05-23
# o Extracted platform-independent AromaTransform from Transform.
# o All Affymetrix-specific code is now in Transform.AFFX.R.
# o Removed some dependencies to CDFs.
# 2007-12-08
# o getOutputDataSet() of Transform was updated to utilize the new 'cdf'
#   argument in static fromFiles() of AffymetrixCelSet.  This way the 
#   default is not queried (in case it does not exist).
# 2007-09-18
# o Now getOutputDataSet() of Transform carry down certain arguments from
#   the input data set. This will speed up things.
# 2007-09-12
# o Now getOutputDataSet() of Transform passes down '...'static
#   fromFiles() of the AffymetrixCelSet class being setup.
# o Now isDone() of Transform throws an error if too many output files are
#   found.  Before it used to return FALSE.
# 2007-09-05
# o Added test against generating an output path that is the same as the
#   path of the input data set.
# 2007-06-25
# o BUG FIX: When getOutputDataSet() retrieved the output data set, the chip
#   type of the CEL files would be validated against the path name, also when
#   then CDF of the input set was overriden.  Now the output data set is
#   setup using 'checkChipType=FALSE'.  Thanks Mark Robinson for 
#   troubleshooting this.
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-02-28
# o Now getOutputData() of Transform make sure to pass down the CDF too.
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2007-01-07
# o BUG FIX: getOutputFiles() would return "private" (prefix '.') files too.
#   This caused for instance FragmentLengthNormalization to return FALSE
#   for isDone() after being average, because one too many files was found.
# 2007-01-06
# o Renamed to Transform (from Preprocessing).
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
