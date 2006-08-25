###########################################################################/**
# @RdocClass AffymetrixDataSet
#
# @title "The AffymetrixDataSet class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixDataSet object represents an Affymetrix data set consisting
#  of one or several Affymetrix data files with \emph{identical} chip types.

#  The purpose of this class is to keep memory usage down.  
#  For this reason, each time data is requested, internally it is read from
#  file and optionally transformed via, say, normalization functions.
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataFiles}{A @list of @see "AffymetrixDataFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @examples "../incl/AffymetrixDataSet.Rex"
#
# \seealso{
#   @see "AffymetrixApdFile" and @see "AffymetrixCelFile".
# }
#
# @author
#*/###########################################################################
setConstructorS3("AffymetrixDataSet", function(dataFiles=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'dataFiles':
  if (is.null(dataFiles)) {
  } else if (is.list(dataFiles)) {
  } else if (inherits(dataFiles, "AffymetrixDataSet")) {
    return(as.AffymetrixDataSet(dataFiles));
  } else {
    throw("Argument 'dataFiles' is of unknown type: ", mode(dataFiles));
  }

  extend(Object(), "AffymetrixDataSet",
    .apdMap = NULL,
    dataFiles = as.list(dataFiles)
  )
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the Affymetrix data set"
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
setMethodS3("as.character", "AffymetrixDataSet", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Path: ", getPath(this), ".", sep="");
  s <- paste(s, " Number of arrays: ", nbrOfArrays(this), ".", sep="");
  s <- paste(s, " Chip type: ", getChipType(this), ".", sep="");
  s <- paste(s, " File type: ", getFileType(this), ".", sep="");
  if (hasTransforms(this)) {
    s <- paste(s, " Has tranforms: TRUE.", sep="");
    s <- paste(s, " Tranform signals: ", isTransforming(this), ".", sep="");
  }
  if (!is.null(getApdMap(this))) {
    s <- paste(s, " Using reading map: ", getApdMap(this), sep="");
  }
  s;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) of the data set"
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
setMethodS3("getPath", "AffymetrixDataSet", function(this, ...) {
  getPath(this$dataFiles[[1]]);
})


###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays in the data set"
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
setMethodS3("nbrOfArrays", "AffymetrixDataSet", function(this, ...) {
  length(this$dataFiles);
})


###########################################################################/**
# @RdocMethod getNames
#
# @title "Gets the names of the arrays in the data set"
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
setMethodS3("getNames", "AffymetrixDataSet", function(this, ...) {
  unlist(lapply(this, FUN=getName))
})




###########################################################################/**
# @RdocMethod getProbesetIntensities
#
# @title "Gets probe signals for a subset of probesets and a subset of arrays"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{arrays}{An @integer index @vector specifying arrays to be read.
#    If @NULL, all arrays are read.}
#  \item{probesets}{An @integer index @vector specifying probesets to be read. 
#    If @NULL, all probesets are read.}
#  \item{readMap}{An @integer @vector specifying a read map used to read the
#    probe signals.}
#  \item{doTransform}{If @TRUE, probe signals are transformed by the
#    internal transformation functions, otherwise not.}
#  \item{...}{Arguments passed to the low-level function for read probesets, 
#    e.g. @see "affxparser::readCelUnits" or @see "aroma.apd::readApdUnits".}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a named @list structure.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbesetIntensities", "AffymetrixDataSet", function(this, probesets=NULL, readMap=getReadMap(this), doTransform=isTransforming(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'doTransform':
  doTransform <- Arguments$getLogical(doTransform);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Transform signals?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfArrays <- nbrOfArrays(this);
  if (doTransform && hasTransforms(this)) {
     # Building transform functions
     verbose && enter(verbose, "Building transform functions");
     transforms <- list("vector", nbrOfArrays);
     for (kk in seq(length=nbrOfArrays)) {
       transforms[[kk]] <- function(y, ...) {
         transformProbeSignals(this, array=kk, y);
       }
     }
     verbose && exit(verbose);
  } else {
     transforms <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pathnames <- unlist(lapply(this, getPathname), use.names=FALSE);

  # Get first data file
  dataFile <- this$dataFiles[[1]];

  # Get the type
  type <- getFileType(this);
  fcnName <- paste("read", capitalize(type), "Units", sep="");
  if (!exists(fcnName, mode="function"))
    throw("Cannot read data by probesets.  Unknown file type: ", type);

  fcn <- get(fcnName, mode="function");
  fcn(pathnames, units=probesets, readMap=readMap, ..., dropArrayDim=FALSE, transforms=transforms, verbose=verbose);
}) # getProbesetIntensities()



setMethodS3("getFileType", "AffymetrixDataSet", function(this, ...) {
  getFileType(this$dataFiles[[1]], ...);
})

setMethodS3("getCdf", "AffymetrixDataSet", function(this, ...) {
  getCdf(this$dataFiles[[1]], ...);
})


###########################################################################/**
# @RdocMethod seq
#
# @title "Gets an vector of data file indices"
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
setMethodS3("seq", "AffymetrixDataSet", function(this, ...) {
  seq(length=nbrOfArrays(this));
})



###########################################################################/**
# @RdocMethod lapply
#
# @title "Applies a function to each of the Affymetrix data files"
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
setMethodS3("lapply", "AffymetrixDataSet", function(this, ...) {
  lapply(this$dataFiles, ...);
})


###########################################################################/**
# @RdocMethod as.list
#
# @title "Returns the data files of the Affymetrix data set"
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
#  Returns a @list of @see "AffymetrixDataFile"s.
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
setMethodS3("as.list", "AffymetrixDataSet", function(x, ...) {
  # To please R CMD check.
  this <- x;

  this$dataFiles;
})



###########################################################################/**
# @RdocMethod extract
#
# @title "Extract a subset of the data set"
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
#   Returns an @see "AffymetrixDataSet" (or a subclass) object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extract", "AffymetrixDataSet", function(this, arrays, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  arrays <- Arguments$getIndices(arrays, range=range(seq(this)));

  res <- clone(this);
  res$dataFiles <- this$dataFiles[arrays];

  res;
})



###########################################################################/**
# @RdocMethod as.AffymetrixDataSet
# @alias as.AffymetrixDataSet.list
# @alias as.AffymetrixDataSet.default
#
# @title "Coerce an object to an AffymetrixDataSet object"
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
#   Returns an @see "AffymetrixDataSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixDataSet", "AffymetrixDataSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixDataSet", "list", function(object, ...) {
  for (kk in seq(along=object)) {
    if (!inherits(object[[kk]], "AffymetrixDataFile")) { 
      throw("List element #", kk, " of argument 'object' is not an AffymetrixDataFile: ", class(object[[kk]])[1]);
    }
  }
  AffymetrixDataSet(object, ...);
})

setMethodS3("as.AffymetrixDataSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixDataSet object: ", mode(object));
})



###########################################################################/**
# @RdocMethod getProbeIntensities
#
# @title "Gets probe intensities from a set of probes and a set of arrays"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{arrays}{An @integer index @vector specifying arrays to be read.
#    If @NULL, all arrays are read.}
#  \item{probes}{An @integer index @vector specifying probes to be read. 
#    If @NULL, all probes are read.}
#  \item{doTransform}{If @TRUE, probe signals are transformed by the
#    internal transformation functions, otherwise not.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric \eqn{NxK} matrix, where \eqn{N} is the number of 
#   probes read, and \eqn{K} is the number of arrays read.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbeIntensities", "AffymetrixDataSet", function(this, arrays=NULL, probes=NULL, doTransform=isTransforming(this), verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probes':
  nbrOfProbes <- nbrOfProbes(this)[1];
  if (!is.null(probes)) {
    probes <- Arguments$getIndices(probes, range=c(1,nbrOfProbes));
    nbrOfProbes <- length(probes);
  }

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq(this);
  } else {
    arrays <- Arguments$getIndices(arrays, range=c(1,nbrOfArrays(this)));
  }

  # Argument 'doTransform':
  doTransform <- Arguments$getLogical(doTransform);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading probe signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(arrays);

  verbose && enter(verbose, "Getting probe signals for ", nbrOfArrays, " arrays.");

  # Allocating the return matrix
  res <- matrix(NA, nrow=nbrOfProbes, ncol=nbrOfArrays);

  for (kk in seq(along=arrays)) {
    verbose && enter(verbose, "Array #", kk, " of ", nbrOfArrays);
    array <- arrays[kk];
    dataFile <- this$dataFiles[[array]];

    res[,kk] <- getProbeIntensities(dataFile, probes=probes, verbose=verbose);

    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  res;
}) # getProbeIntensities()



###########################################################################/**
# @RdocMethod "[["
#
# @title "Gets a subset of the data files in the data set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{i}{An @integer index @vector specifying data files to be returned.}
#  \item{drop}{If @TRUE and only one array is returned, the data file is
#    return directly, otherwise as an element in a list.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list of @see "AffymetrixDataFile" objects, or a single
#   @see "AffymetrixDataFile" object if \code{drop} is @TRUE and only
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
setMethodS3("[[", "AffymetrixDataSet", function(this, i=NULL, drop=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'i':
  if (is.null(i))
    return(this$dataFiles);
  i <- Arguments$getIndices(i, range=range(seq(this)));

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);


  files <- this$dataFiles[i];
  if (length(i) == 1 && drop)
    files <- files[[1]];
  files;
})


###########################################################################/**
# @RdocMethod "["
#
# @title "Extract a subset of the data set"
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
setMethodS3("[", "AffymetrixDataSet", function(this, i=NULL, ...) {
  extract(this, arrays=i, ...);
})


setMethodS3("clearCache", "AffymetrixDataSet", function(this, ...) {
  # Clear the cache of all data files
  lapply(this, clearCache);

  # Clear the cache of the CDF object
  clearCache(getCdf(this));

  # Then for this object
  NextMethod("clearCache", this);
})


# setMethodS3("[", "AffymetrixDataSet", function(this, i=NULL, j=NULL, drop=FALSE) {
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   # Validate arguments
#   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   # Argument 'i':
#   if (is.null(i))
#     i <- NULL;
#   # Argument 'j':
#   if (is.null(j))
#     j <- NULL;
# 
#   y <- getProbeIntensities(this, probes=i, arrays=j);
#   if (drop && length(j) == 1) {
#     y <- drop(y);
#   }
# 
#   y;
# })
# 




############################################################################
# HISTORY:
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
# o Added argument 'subset' to calcAvgProbeSignals() & normalizeQuantile().
# 2006-03-15
# o Now nbrOfProbes() returns the number of probes for the first file only.
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
