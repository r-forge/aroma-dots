###########################################################################/**
# @RdocClass AffymetrixFileSet
#
# @title "The AffymetrixFileSet class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixFileSet object represents a set of @see "AffymetrixFile"s
#  with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AffymetrixFileSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    lapply(files, FUN=function(df) {
      if (!inherits(df, "AffymetrixFile"))
        throw("Argument 'files' contains a non-AffymetrixFile object: ", class(df));
    })
  } else if (inherits(files, "AffymetrixFileSet")) {
    return(as.AffymetrixFileSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }

  extend(Object(), "AffymetrixFileSet",
    files = as.list(files)
  )
})



setMethodS3("getName", "AffymetrixFileSet", function(this, ...) {
  # The name of a file set is inferred from the pathname of the directory
  # of the set assuming path/to/<name>/<something>/<"chip type">/
  # Get the path of this file set
  path <- getPath(this);

  # path/to/<name>/<something>
  path <- dirname(path);

  # path/to/<name>
  path <- dirname(path);

  # <name>
  name <- basename(path);

  name;
})


###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the Affymetrix file set"
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
setMethodS3("as.character", "AffymetrixFileSet", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Number of files: %d", nbrOfFiles(this)));
  s <- c(s, sprintf("Total file size: %.2fMb", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})

setMethodS3("getFileSize", "AffymetrixFileSet", function(this, ...) {
  if (is.null(fileSize <- this$.fileSize)) {
    fileSize <- sum(unlist(lapply(this, FUN=getFileSize), use.names=FALSE));
    this$.fileSize <-  fileSize;
  }
  fileSize;
})

###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) of the file set"
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
setMethodS3("getPath", "AffymetrixFileSet", function(this, ...) {
  getPath(this$files[[1]]);
})



###########################################################################/**
# @RdocMethod length
# @aliasmethod nbrOfFiles
#
# @title "Gets the number of files in the file set"
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
setMethodS3("length", "AffymetrixFileSet", function(this, ...) {
  length(this$files);
})

setMethodS3("nbrOfFiles", "AffymetrixFileSet", function(this, ...) {
  length(this, ...);
})



###########################################################################/**
# @RdocMethod lapply
#
# @title "Applies a function to each of the Affymetrix file files"
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
setMethodS3("lapply", "AffymetrixFileSet", function(this, ...) {
  lapply(this$files, ...);
})


###########################################################################/**
# @RdocMethod getNames
#
# @title "Gets the names of the files in the file set"
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
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getNames", "AffymetrixFileSet", function(this, ...) {
  unlist(lapply(this, FUN=getName))
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
#   Returns an @integer @vector in [1,N] where N is the number of files,
#   or an empty vector if the file set is empty.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("seq", "AffymetrixFileSet", function(this, ...) {
  seq(length=nbrOfFiles(this));
})


###########################################################################/**
# @RdocMethod as.list
#
# @title "Returns the files of the Affymetrix file set"
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
#  Returns a @list of @see "AffymetrixFile"s.
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
setMethodS3("as.list", "AffymetrixFileSet", function(x, ...) {
  # To please R CMD check.
  this <- x;

  this$files;
})



###########################################################################/**
# @RdocMethod extract
#
# @title "Extract a subset of the file set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{files}{File indices.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" (or a subclass) object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extract", "AffymetrixFileSet", function(this, files, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  files <- Arguments$getIndices(files, range=range(seq(this)));

  res <- clone(this);
  res$files <- this$files[files];
  clearCache(res);  # Some cached values are incorrect now.

  res;
})


###########################################################################/**
# @RdocMethod as.AffymetrixFileSet
# @alias as.AffymetrixFileSet.list
# @alias as.AffymetrixFileSet.default
#
# @title "Coerce an object to an AffymetrixFileSet object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixFileSet", "AffymetrixFileSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixFileSet", "list", function(object, ...) {
  AffymetrixFileSet(object, ...);
})

setMethodS3("as.AffymetrixFileSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixFileSet object: ", mode(object));
})


setMethodS3("clearCache", "AffymetrixFileSet", function(this, ...) {
  # Clear the cache of all files
  lapply(this, clearCache);

  # Clear the cache of the CDF object
  clearCache(getCdf(this));

  # Then for this object
  NextMethod("clearCache", this);
})




###########################################################################/**
# @RdocMethod fromFiles
#
# @title "Defines an AffymetrixFileSet object by searching for Affymetrix files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{path}{The path where to search for Affymetrix files.}
#  \item{pattern}{The filename pattern for match files. 
#     If @NULL, filename extensions corresponding to known subclasses
#     of the abstract @see "AffymetrixFile" class are search for.}
#  \item{recursive}{If @TRUE, subdirectories are search recursively,
#     otherwise not.}
#  \item{...}{Other arguments passed to @see "base::list.files"}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" object.
# }
#
# \section{Reserved filenames}{
#   Note that files with names starting with a period \code{.} are not 
#   searched for.  The reason for this is that such files are reserved for
#   internal use of this package.  For instance, the package store average
#   signals across CEL files in a file named as \code{.average<something>.CEL}
#   in the same directory as the CEL files of the dataset, and when such a
#   directory is scanned we do not want such files to be interpreted as data.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromFiles", "AffymetrixFileSet", function(static, path=NULL, pattern=NULL, recursive=TRUE, fileClass="AffymetrixFile", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'pattern':
  if (!is.null(pattern))
    pattern <- Arguments$getRegularExpression(pattern);

  # Argument 'fileClass':
  clazz <- Class$forName(fileClass);
  dfStatic <- getStaticInstance(clazz);
  if (!inherits(dfStatic, "AffymetrixFile"))
    throw("Argument 'fileClass' is not refering to a AffymetrixFile class: ", fileClass);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create AffymetrixFile object from the matching files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan for Affymetrix files
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE, 
                               all.files=FALSE, recursive=recursive, ...);
  if (length(pathnames) == 0)
    throw("No files found: ", path);

  files <- list();
  for (kk in seq(along=pathnames)) {
    df <- fromFile(dfStatic, pathnames[kk]);
    files[[kk]] <- df;
    if (kk == 1) {
      clazz <- Class$forName(class(df)[1]);
      dfStatic <- getStaticInstance(clazz);
    }
  }

  # Create the file set object
  set <- newInstance(static, files);

  set;
}, static=TRUE)



############################################################################
# HISTORY:
# 2006-08-27
# o Made filenames starting with a period reserved for internal use.
# 2006-08-26
# o Now getName() of a file set is inferred from the pathname:
#     path/to/<name>/chip_files/<"chip type">/
# 2006-08-21
# o Renamed 'array' to 'file'.
# o Extracted from AffymetrixCelSet.R.
# 2006-08-11
# o Added clearCache() which also clears the cache of all data file object.
# 2006-05-16
# o Redefined "[" to extract arrays.
# 2006-04-13
# o Added Rdoc comments for all methods.
# 2006-04-09
# o Now the read map is loaded automatically when fromFiles() used.
# 2006-03-30
# o Updated to new aroma.apd.
# 2006-03-18
# o Added argument 'subset' to calcAvgCellSignals() & normalizeQuantile().
# 2006-03-15
# o Now nbrOfCells() returns the number of cells for the first file only.
# o Now the fromFiles(static, ...) creates an object of the same class as 
#   the static object.
# 2006-03-04
# o Added mapping functions.
# o Added writeApd().
# 2006-03-03
# o Added lapply().
# 2006-03-02
# o Updated to deal with AffymetrixDataFile object instead of CEL directly.
# 2006-02-21
# o Letting readCelUnits() transform signals improves speed substantially.
# o Making use of new multi-array readCelUnits().
# 2006-02-20
# o Created.
############################################################################
