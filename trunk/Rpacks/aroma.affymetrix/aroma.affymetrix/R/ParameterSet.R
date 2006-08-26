###########################################################################/**
# @RdocClass ParameterSet
#
# @title "The abstract ParameterSet class"
#
# \description{
#  @classhierarchy
#
#  An ParameterSet object represents a set of parameter estimates.
# }
# 
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the data file.}
#   \item{path}{An optional path to the data data file.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
# @visibility "private"
#*/###########################################################################
setConstructorS3("ParameterSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }

  extend(Object(), "ParameterSet",
    files = files
  )
}, private=TRUE)



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing this parameter set"
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
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "ParameterSet", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Number of arrays: ", nbrOfArrays(this), ".", sep="");
  s <- paste(s, " Name: ", getName(this), ".", sep="");
  s;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) of this parameter set"
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
setMethodS3("getPath", "ParameterSet", function(this, ...) {
  getPath(this$files[[1]]);
})




###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of this parameter set"
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
setMethodS3("getName", "ParameterSet", function(this, ...) {
  basename(getPath(this));
})



###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays in this parameter set"
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
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfArrays", "ParameterSet", function(this, ...) {
  length(this$files);
})





###########################################################################/**
# @RdocMethod seq
#
# @title "Gets an vector of file indices"
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
#   Returns an @integer @vector in [1,N] where N is the number of arrays,
#   or an empty vector if the data set is empty.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("seq", "ParameterSet", function(this, ...) {
  seq(length=nbrOfArrays(this));
})



###########################################################################/**
# @RdocMethod lapply
#
# @title "Applies a function to each of the parameter files"
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
#   Returns a @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("lapply", "ParameterSet", function(this, ...) {
  lapply(this$files, ...);
})



###########################################################################/**
# @RdocMethod as.list
#
# @title "Returns the data files of the parameter set"
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
#  Returns a @list of @see "ParameterFile"s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.list", "ParameterSet", function(x, ...) {
  # To please R CMD check.
  this <- x;

  this$files;
})



###########################################################################/**
# @RdocMethod extract
#
# @title "Extract a subset of the parameter set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{arrays}{Array indices.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "ParameterSet" (or a subclass) object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extract", "ParameterSet", function(this, arrays, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  arrays <- Arguments$getIndices(arrays, range=range(seq(this)));

  res <- clone(this);
  res$files <- this$files[arrays];

  res;
})


###########################################################################/**
# @RdocMethod "[["
#
# @title "Gets a subset of the files in the parameter set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{i}{An @integer index @vector specifying files to be returned.}
#  \item{drop}{If @TRUE and only one array is returned, the data file is
#    return directly, otherwise as an element in a list.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list of @see "ParameterFile" objects, or a single
#   @see "ParameterFile" object if \code{drop} is @TRUE and only
#   one data file was selected.
# }
#
# @author
#
# \seealso{
#   @seemethod "[".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("[[", "ParameterSet", function(this, i=NULL, drop=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'i':
  if (is.null(i))
    return(this$files);
  i <- Arguments$getIndices(i, range=range(seq(this)));

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);


  files <- this$files[i];
  if (length(i) == 1 && drop)
    files <- files[[1]];
  files;
})


###########################################################################/**
# @RdocMethod "["
#
# @title "Extract a subset of the parameter set"
#
# \description{
#   @get "title".
#   This is just a wrapper for @seemethod "extract".
# }
#
# @synopsis
#
# \arguments{
#  \item{i}{An @integer index @vector specifying the arrays to be returned.}
#  \item{...}{Arguments passed to @seemethod "extract".}
# }
#
# \value{
#   Returns what @seemethod "extract" returns.
# }
#
# @author
#
# \seealso{
#   @seemethod "extract".
#   @seemethod "[[".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("[", "ParameterSet", function(this, i=NULL, ...) {
  extract(this, arrays=i, ...);
})



###########################################################################/**
# @RdocMethod fromFiles
#
# @title "Defines a ParameterSet object by searching for files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{path}{The path where to search for files.}
#  \item{pattern}{The filename pattern for match files. 
#     If @NULL, filename extensions corresponding to known subclasses
#     of the abstract @see "ParameterFile" class are search for.}
#  \item{recursive}{If @TRUE, subdirectories are search recursively,
#     otherwise not.}
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "ParameterSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromFiles", "ParameterSet", function(static, path="estimates/", pattern=NULL, recursive=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Get the Class object of the ParameterFile class corresponding to
  # the static instance 'static'.
  className <- gsub("Set$", "File", class(static)[1]);
  clazz <- Class$forName(className);

  # Argument 'pattern':
  if (is.null(pattern)) {
    # Build regular expression pattern from known subclasses
    classNames <- c(className, getKnownSubclasses(clazz));

    # File name extensions to search for
    pattern <- gsub("Parameter", "Parameter([A-Z][a-zA-Z0-9]*)", className);
    pattern <- paste("^", pattern, "$", sep="");

    exts <- c();
    for (name in classNames) {
      clazz <- Class$forName(name);
      ext <- clazz$getExtension();
      exts <- c(exts, ext);
    }
    exts <- unique(exts);
    pattern <- paste("[.](", paste(exts, collapse="|"), ")$", sep="");
  } else {
    pattern <- Arguments$getRegularExpression(pattern);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create ParameterFile object from the matching files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan for Affymetrix data files files
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE, 
                                                  recursive=recursive, ...);
  if (length(pathnames) == 0)
    throw("No data files found: ", path);
  files <- list();
  dfStatic <- getStaticInstance(clazz);
  for (kk in seq(along=pathnames)) {
    df <- fromFile(dfStatic, pathnames[kk]);
    files[[kk]] <- df;
    if (kk == 1) {
      clazz <- Class$forName(class(df)[1]);
      dfStatic <- getStaticInstance(clazz);
    }
  }

  # Create the data-set object
  newInstance(static, files);
}, static=TRUE)






############################################################################
# HISTORY:
# 2006-05-31
# o Created.  Eventually there might be superclasses AffymetrixFileSet and
#   AffymetrixFile.
############################################################################
