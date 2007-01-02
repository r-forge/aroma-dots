###########################################################################/**
# @RdocClass AffymetrixCelSet
#
# @title "The AffymetrixCelSet class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCelSet object represents a data set of Affymetrix CEL files 
#  with \emph{identical} chip types.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixCelFile":s.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \examples{\dontrun{
#   @include "../incl/AffymetrixCelSet.Rex"
# }}
#
# \seealso{
#   @see "AffymetrixCelFile".
# }
#
# @author
#*/###########################################################################
setConstructorS3("AffymetrixCelSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    lapply(files, FUN=function(df) {
      if (!inherits(df, "AffymetrixCelFile"))
        throw("Argument 'files' contains a non-AffymetrixCelFile object: ", class(df));
    })
  } else if (inherits(files, "AffymetrixCelSet")) {
    return(as.AffymetrixCelSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }


  extend(AffymetrixFileSet(files=files, ...), "AffymetrixCelSet",
    "cached:.intensities" = NULL,
    "cached:.intensitiesIdxs" = NULL,
    "cached:.unitsCache" = NULL,
    "cached:.getUnitIntensitiesCache" = NULL,
    "cached:.averageFiles" = list(),
    "cached:.timestamps" = NULL,
    "cached:.fileSize" = NULL
  )
})


setMethodS3("clone", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Cloning Affymetrix CEL set");

  # Clone itself and the files.  The call below will clear the cache!
  object <- NextMethod("clone", clear=FALSE, ..., verbose=less(verbose));

  # Clone the CDF (this will update the CDF of all file object)
  verbose && enter(verbose, "Cloning CDF");
  cdf <- clone(getCdf(object));
  verbose && exit(verbose);
  verbose && enter(verbose, "Adding CDF to CEL set");
  setCdf(object, cdf, .checkArgs=FALSE);
  verbose && exit(verbose);

  verbose && exit(verbose);

  object;
})


setMethodS3("append", "AffymetrixCelSet", function(this, other, clone=TRUE, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (!inherits(other, class(this)[1])) {
    throw("Argument 'other' is not an ", class(this)[1], " object: ", 
                                                      class(other)[1]);
  }

  verbose && enter(verbose, "Appending CEL set");
  verbose && print(verbose, other);

  # Validate chip type
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  for (file in getFiles(other)) {
    oCdf <- getCdf(file);
    oChipType <- getChipType(oCdf);
    if (!identical(oChipType, chipType)) {
      throw("Argument 'other' contains a CEL file of different chip type: ",
                                                oChipType, " != ", chipType);
    }
  }

  # Append other
  this <- NextMethod("append", this, other=other, clone=clone, ...);

  # Set the same CDF for all CEL files
  verbose && enter(verbose, "Updating the CDF for all files");
  setCdf(this, cdf);
  verbose && exit(verbose);

  verbose && exit(verbose);

  this;
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
setMethodS3("as.character", "AffymetrixCelSet", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  s <- c(s, sprintf("Tags: %s", tags));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("Number of arrays: %d", nbrOfArrays(this)));
  # Get CEL header timestamps
  ts <- getTimestamps(this);
  ts <- range(ts);
  ts <- format(ts, "%Y-%m-%d %H:%M:%S");  # range() gives strange values?!?
  s <- c(s, sprintf("Time period: %s -- %s", ts[1], ts[2]));
  s <- c(s, sprintf("Total file size: %.2fMb", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("getTimestamps", "AffymetrixCelSet", function(this, ..., force=FALSE) {
  ts <- this$.timestamps;

  if (force || is.null(ts)) {
    # Get CEL header dates
    ts <- lapply(this, getTimestamp);
    ts <- do.call("c", args=ts);
    this$.timestamps <- ts;
  }

  ts;
})

setMethodS3("getIdentifier", "AffymetrixCelSet", function(this, ..., force=FALSE) {
  identifier <- this$.identifier;
  if (force || is.null(identifier)) {
    identifier <- NextMethod("getIdentifier");
    if (is.null(identifier)) {
      identifiers <- lapply(this, getIdentifier);
      identifier <- digest(identifiers);
    }
    this$.identifier <- identifier;
  }
  identifier;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getSampleNames
#
# @title "Gets the names of the samples in the file set"
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
setMethodS3("getSampleNames", "AffymetrixCelSet", function(this, ...) {
  unlist(lapply(this, FUN=getSampleName))
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getSiblings
#
# @title "Gets the all data sets that refers to the same samples a this one"
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
#  Returns a named @list of @see "AffymetrixCelSet" objects.
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
setMethodS3("getSiblings", "AffymetrixCelSet", function(this, notSelf=FALSE, ...) {
  # Scan parent directory for all possible data sets.
  # /path/to/<data set name>/chip_data/<chip type>
  path <- getPath(this);

  # /path/to/<data set name>/chip_data
  parent <- getParent(path);
  # /path/to/<data set name>
  parent <- getParent(parent);

  # /path/to/
  dataPath <- getParent(parent);

  # Now scan this file directory tree
  paths <- findCelSet(name=getName(this), paths=dataPath, firstOnly=FALSE);

  if (notSelf)
    paths <- paths[(paths != getPath(this))];

  sets <- vector("list", length(paths));
  names(sets) <- basename(paths);
  for (kk in seq(along=paths)) {
    path <- paths[kk];
    if (!notSelf && identical(path, getPath(this))) {
      sets[[kk]] <- this;
    } else {
      sets[[kk]] <- fromFiles(this, path=path);
    }
  }
  
  sets;
})


###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CEL set"
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
setMethodS3("getCdf", "AffymetrixCelSet", function(this, ...) {
  getCdf(this$files[[1]], ...);
})


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this CEL set"
#
# \description{
#  @get "title".  This structures is assigned to all CEL files in the set.
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setCdf", "AffymetrixCelSet", function(this, cdf, verbose=FALSE, ...) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Setting CDF for data set");
  verbose && print(verbose, cdf);

  # Nothing to do?
  oldCdf <- getCdf(this);
#  if (equals(cdf, oldCdf))
#    return(invisible(this));

  # Set the CDF for all CEL files
  verbose && enter(verbose, "Setting CDF for each data file");
  lapply(this, setCdf, cdf, ...);
  verbose && exit(verbose);

  # Have to clear the cache 
  verbose && enter(verbose, "Clearing data-set cache");
  clearCache(this);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(this);
})



setMethodS3("fromFiles", "AffymetrixCelSet", function(static, path="chip_data/", pattern="[.](c|C)(e|E)(l|L)$", ..., fileClass="AffymetrixCelFile", verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Defining ", class(static)[1], " from files");

  this <- fromFiles.AffymetrixFileSet(static, path=path, pattern=pattern, ..., fileClass=fileClass, verbose=less(verbose));

  # Use the same CDF object for all CEL files.
  verbose && enter(verbose, "Updating the CDF for all files");
  setCdf(this, getCdf(this), .checkArgs=FALSE);
  verbose && exit(verbose);

  verbose && exit(verbose);

  this;
})




###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays in the file set"
#
# \description{
#   @get "title".
#   This is just a wrapper for \code{nbrOfFiles()}.
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
setMethodS3("nbrOfArrays", "AffymetrixCelSet", function(this, ...) {
  nbrOfFiles(this, ...);
})



###########################################################################/**
# @RdocMethod as.AffymetrixCelSet
# @alias as.AffymetrixCelSet.list
# @alias as.AffymetrixCelSet.default
#
# @title "Coerce an object to an AffymetrixCelSet object"
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
#   Returns an @see "AffymetrixCelSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixCelSet", "AffymetrixCelSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixCelSet", "list", function(object, ...) {
  AffymetrixCelSet(object, ...);
})

setMethodS3("as.AffymetrixCelSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixCelSet object: ", mode(object));
})



###########################################################################/**
# @RdocMethod isDuplicated
#
# @title "Identifies duplicated data files"
#
# \description{
#   @get "title" by comparing the timestamps in the CEL headers.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @logical @vector of length equal to the number of files
#   in the data set.
#   An element with value @TRUE indicates that the corresponding data file
#   has the same time stamp as another preceeding data file.
# }
#
# \examples{\dontrun{
#   # The data set of interest
#   ds <- AffymetrixCelSet$fromFiles(path=...)
#
#   # Added other data sets to be used as a reference
#   for (path in refPaths) {
#     dsR <- AffymetrixCelSet$fromFiles(path=path)
#     append(ds, dsR)
#   }
#
#   # Keep only unique arrays
#   ds <- extract(ds, !isDuplicated(ds))
# }}
#
# @author
#
# \seealso{
#   Internally @see "base::duplicated" is used to compare timestamps.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDuplicated", "AffymetrixCelSet", function(this, ...) {
  # Get the CEL header timestamp for all files
  timestamps <- getTimestamps(this);

  dups <- duplicated(timestamps);
  names(dups) <- getNames(this);

  dups;
})



setMethodS3("getData", "AffymetrixCelSet", function(this, indices=NULL, fields=c("x", "y", "intensities", "stdvs", "pixels"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  nbrOfCells <- nbrOfCells(getCdf(this));
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices, range=c(1, nbrOfCells));
    nbrOfCells <- length(indices);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  nbrOfArrays <- nbrOfArrays(this);
  verbose && enter(verbose, "Getting cell data for ", nbrOfArrays, " arrays.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating the return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating the return structure");
  nbrOfFields <- length(fields);
  res <- vector("list", nbrOfFields);
  names(res) <- fields;
  for (field in fields) {
    res[[field]] <- matrix(NA, nrow=nbrOfCells, ncol=nbrOfArrays);
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading cell signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Requiring data from ", nbrOfArrays, " arrays");
  for (kk in seq(length=nbrOfArrays)) {
    verbose && enter(verbose, "Array #", kk, " of ", nbrOfArrays);
    dataFile <- this$files[[kk]];
    value <- getData(dataFile, indices=indices, fields=fields, verbose=less(verbose));
    for (field in fields) {
      res[[field]][,kk] <- value[[field]];
      value[[field]] <- NULL;
    }
    rm(value); gc();
    verbose && exit(verbose);
  }
  verbose && exit(verbose);


  verbose && exit(verbose);

  res;
}) # getData()


###########################################################################/**
# @RdocMethod getIntensities
#
# @title "Gets cell intensities from a set of cells and a set of arrays"
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
#  \item{indices}{An @integer index @vector specifying cells to be read. 
#    If @NULL, all cells are read.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric \eqn{NxK} matrix, where \eqn{N} is the number of 
#   cells read, and \eqn{K} is the number of arrays read.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getIntensities", "AffymetrixCelSet", function(this, ...) {
  getData(this, ..., fields="intensities")$intensities;
}) # getIntensities()



###########################################################################/**
# @RdocMethod getUnitIntensities
#
# @title "Gets cell signals for a subset of units and a subset of arrays"
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
#  \item{units}{An @integer index @vector specifying units to be read. 
#    If @NULL, all units are read.}
#  \item{readMap}{An @integer @vector specifying a read map used to read the
#    cell signals.}
#  \item{...}{Arguments passed to the low-level function for read units, 
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
setMethodS3("getUnitIntensities", "AffymetrixCelSet", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- digest(list(units=units, ...));
  res <- this$.getUnitIntensitiesCache[[key]];
  if (!is.null(res)) {
    verbose && cat(verbose, "getUnitIntensitiesCache(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the pathnames of all CEL files
  pathnames <- unlist(lapply(this, getPathname), use.names=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Cached values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.list(units)) {
    res <- readCelUnits(pathnames, cdf=units, readStdvs=FALSE, 
                                                    readPixels=FALSE, ...);
  } else {
    # Always ask for CDF information from the CDF object!
    cdf <- readUnits(getCdf(this), units=units);
    res <- readCelUnits(pathnames, cdf=cdf, readStdvs=FALSE, 
                                                    readPixels=FALSE, ...);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "readUnits(): Updating cache");
  this$.getUnitIntensitiesCache <- list();
  this$.getUnitIntensitiesCache[[key]] <- res;

  res;
})



setMethodS3("readUnits", "AffymetrixCelSet", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "readCelUnits() of AffymetrixCelSet");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating hashcode key for cache");
  if (is.list(units)) {
    key <- digest(list(units=names(units), ...));
  } else {
    key <- digest(list(units=units, ...));
  }
  verbose && exit(verbose);
  if (!force) {
    verbose && enter(verbose, "Trying to obtain cached data");
    res <- this$.readUnitsCache[[key]];
    verbose && exit(verbose);
    if (!is.null(res)) {
      verbose && cat(verbose, "readUnits(): Returning cached data");
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the pathnames of all CEL files
  pathnames <- getPathnames(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read data from file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calling readCelUnits() for ", 
                                              length(pathnames), " files");
  if (is.list(units)) {
    res <- readCelUnits(pathnames, cdf=units, ...);
  } else {
    # Always ask for CDF information from the CDF object!
     verbose && enter(verbose, "Retrieving CDF unit information");
    cdf <- readUnits(getCdf(this), units=units, verbose=less(verbose));
    verbose && exit(verbose);
    res <- readCelUnits(pathnames, cdf=cdf, ...);
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  this$.readUnitsCache <- list()
  this$.readUnitsCache[[key]] <- res
  verbose && cat(verbose, "readUnits(): Updated cache");

  verbose && exit(verbose);

  res;
})


###########################################################################/**
# @RdocMethod getAverageFile
#
# @title "Calculates the mean and the standard deviation cell signal across the data set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{The label of the calculated parameters.
#    If @NULL, a default name format \code{<prefix>-<mean>-<sd>} is used.}
#  \item{indices}{An @integer @vector specifying which cells to consider.
#    If \code{"remaining"}, only parameters for cells that have not been
#    are calculated.
#    If @NULL, all cells are used.}
#  \item{mean}{A @character of a @function specifying the function used
#    to calculate the average.}
#  \item{sd}{A @character of a @function specifying the function used
#    to calculate the standard deviation.}
#  \item{na.rm}{If @TRUE, @NAs are excluded before, otherwise not.}
#  \item{...}{Not used.}
#  \item{cellsPerChunk}{A @integer specifying the total number of cells 
#    (across arrays) read into memory per chunk.}
#  \item{moreCells}{A @double scalar indicating if more or less cells
#    should be used per chunk.}
#  \item{force}{If @TRUE, parameters for cells already calculated are
#    recalculated, otherwise not.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns an @see "AffymetrixCelSet" of the same class as the data set
#   averaged.
# }
#
# \details{
#   The parameter estimates are stored as a CEL file of the same class as
#   the data files in the data set.  The CEL file is named \code{<name>.cel}
#   and placed in the directory of the data set.
#   Currently there is no specific data class for this file, but the average
#   cell signals are stored as "intensities", the standard deviation of the
#   cell signals as "stddevs", and the number of data points used for each
#   estimate is stored as "pixels".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAverageFile", "AffymetrixCelSet", function(this, name=NULL, prefix="average", indices="remaining", mean=c("median", "mean"), sd=c("mad", "sd"), na.rm=FALSE, g=NULL, h=NULL, ..., cellsPerChunk=moreCells*10^7/length(this), moreCells=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ("median" %in% mean || "mad" %in% sd) {
    # rowMedians():
    if (require(R.native)) {
      rowMedians <- R.native::rowMedians;
    } else {
      # About 3-10 times slower than rowMedians()
      rowMedians <- function(X, ...) {
        apply(X, MARGIN=1, FUN=median, ...);
      }
    }

    # rowMads():
    rowMads <- function(X, centers=rowMedians(X, ...), constant=1.4826, ...) {
      constant * rowMedians(abs(X - centers), ...);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mean':
  if (is.character(mean)) {
    mean <- match.arg(mean);
    meanName <- mean;
    if (mean == "mean") {
      mean <- base::rowMeans;
    } else if (mean == "median") {
      mean <- rowMedians;
    }
  } else if (is.function(mean)) {
    meanName <- "customMean";
  } else {
    throw("Argument 'mean' must be either a character or a function: ", mode(mean));
  }

  # Argument 'sd':
  if (is.character(sd)) {
    sd <- match.arg(sd);
    sdName <- sd;
    if (sd == "sd") {
      sd <- rowSds;
    } else if (sd == "mad") {
      sd <- rowMads;
    }
  } else if (is.function(sd)) {
    sdName <- "customSd";
  } else {
    throw("Argument 'sd' must be either a character or a function: ", 
                                                           mode(sd));
  }

  # Argument 'name':
  if (is.null(name)) {
    key <- list(arrays=sort(getNames(this)), meanName=meanName, sdName=sdName);
    id <- digest(key);
    name <- sprintf("%s-%s-%s,%s", prefix, meanName, sdName, id);
  }

  # Argument 'indices':
  df <- as.list(this)[[1]];
  nbrOfCells <- getHeader(df)$total;
  if (force) {
    if (identical(indices, "remaining")) {
      indices <- NULL;
    }
  }

  if (is.null(indices)) {
    indices <- 1:nbrOfCells; 
  } else if (identical(indices, "remaining")) {
  } else {
    indices <- Arguments$getIndices(indices, range=c(1, nbrOfCells));
  }

  # Argument 'cellsPerChunk':
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving average cell intensities across ", length(this), " arrays");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create CEL file to store the average array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a private filename (with a dot prefix) to make sure it is not
  # identified as a regular CEL file when the directory is scanned for files.
  filename <- sprintf(".%s.CEL", name);
  res <- this$.averageFiles[[filename]];
  if (is.null(res)) {
    res <- createFrom(df, filename=filename, path=getPath(this), verbose=less(verbose));
    this$.averageFiles[[filename]] <- res;
  }

  pathname <- getPathname(res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify which indices to use
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(indices, "remaining")) {
    pixels <- readCel(pathname, readIntensities=FALSE, readStdvs=FALSE, 
                      readPixels=TRUE)$pixels;
    indices <- which(pixels == 0);
    rm(pixels); # Not needed anymore.
  }

  nbrOfIndices <- length(indices);

  # Nothing more to do?
  if (nbrOfIndices == 0)
    return(res);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate the mean and standard deviation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Since we might want to do this robustly, but also because we want to
  # estimate the standard deviation, for each cell we need all data across 
  # arrays at once.  In order to this efficiently, we do this in chunks
  idxs <- 1:nbrOfIndices;
  head <- 1:cellsPerChunk;
  nbrOfChunks <- ceiling(nbrOfIndices / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  # Get the pathnames of all CEL files to average
  pathnames <- lapply(this, getPathname);
  pathnames <- unlist(pathnames, use.names=FALSE);

  if (!na.rm)
    n <- rep(length(this), length=cellsPerChunk);
  count <- 1;
  while (length(idxs) > 0) {
    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < cellsPerChunk) {
      head <- 1:length(idxs);
      if (!na.rm)
        n <- rep(length(pathnames), length=length(idxs));
    }

    # The indices to be used in this chunk
    ii <- idxs[head];
    verbose && cat(verbose, "Chunk size: ", length(ii));

    verbose && enter(verbose, "Reading data");
    X <- readCelIntensities(pathnames, indices=indices[ii]);
    verbose && exit(verbose);

    if (!is.null(g)) {
      verbose && enter(verbose, "Transforming data using y = g(x)");
      X <- g(X);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Estimating averages and standard deviations");
    if (na.rm)
      n <- apply(X, MARGIN=1, FUN=function(x) { sum(!is.na(x)) });

    # Calculate the mean signal    
    mu <- mean(X, na.rm=na.rm);          # Special mean()!
    # Calculate the standard deviation of the signals
    sigma <- sd(X, mean=mu, na.rm=na.rm);   # Special sd()!

    verbose && exit(verbose);

    if (!is.null(h)) {
      verbose && enter(verbose, "Back-transforming estimates using x = h(y)");
      mu <- h(mu);
      sigma <- h(sigma);
      verbose && exit(verbose);
    }

    # Write estimates to result file
    verbose && enter(verbose, "Writing estimates");
    updateCel(pathname, indices=indices[ii], intensities=mu, stdvs=sigma, pixels=n);
    verbose && exit(verbose);

    # Not needed anymore
    mu <- sigma <- NULL;

    # Next chunk...
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    gc();
    verbose && exit(verbose);
  } # while()

  verbose && exit(verbose);

  res;  
})


setMethodS3("getAverage", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, ...);
})

setMethodS3("getAverageLog", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, g=log2, h=function(x) 2^x, ...);
})

setMethodS3("getAverageAsinh", "AffymetrixCelSet", function(this, ...) {
  getAverageFile(this, g=asinh, h=sinh, ...);
})



setMethodS3("range", "AffymetrixCelSet", function(this, ...) {
  range(unlist(lapply(this, FUN=range, ...), use.names=FALSE));
})



setMethodS3("applyToUnitIntensities", "AffymetrixCelSet", function(this, units=NULL, FUN, stratifyBy="pm", verbose=FALSE, ...) {
  y <- getUnitIntensities(this, units=units, stratifyBy=stratifyBy, ...);

  y <- lapply(y, FUN=function(unit) {
    groups <- lapply(unit, FUN=function(group) {
      FUN(group[[1]], ...)
    })
    groups;
  })
  y;
})


setMethodS3("[", "AffymetrixCelSet", function(this, units=NULL, ..., drop=FALSE) {
  res <- readUnits(this, units=units, ...);
  if (drop && length(res) == 1)
    res <- res[[1]];
  res;
})

setMethodS3("[[", "AffymetrixCelSet", function(this, units=NULL, ...) {
  this[units=units, ..., drop=TRUE];
})

###########################################################################/**
# @RdocMethod gcrmaSummary
#
# @title "Calculates the GCRMA expression summary"
#
# \description{
#  @get "title".
#
# Applies GC-based background correction, followed by quantile normalisation
# and then a summary by either median polish or robust regression.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path to save the expression summaries.}
#   \item{name}{Name of the data set containing the expression summaries.}
#   \item{bgPath}{Directory in which to store background-adjusted signals.}
#   \item{normPath}{Directory in which to store normalised signals.}
#   \item{summaryMethod}{Either "medianpolish" or "rlm".}
#   \item{type}{The type of background correction.  Currently accepted types
#       are "fullmodel" (the default, uses MMs) and "affinities" (uses
#       probe sequence only).}
#   \item{indicesNegativeControl}{Locations of any negative control
#       probes (e.g., the anti-genomic controls on the human exon array).
#       If @NULL and type=="affinities", MMs are used as the negative
#       controls.}
#   \item{opticalAdjust}{If @TRUE, apply correction for optical effect,
#       as in @see "gcrma::bg.adjust.optical".}
#   \item{gsbAdjust}{Should we adjust for specific binding (defaults to
#        @TRUE)?}
#   \item{k}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{rho}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{stretch}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{fast}{If @TRUE, an ad hoc transformation of the PM is performed
#       (\code{gcrma::gcrma.bg.transformation.fast}).}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet" containing the expression summaries.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#
# \seealso{
#  @see "gcrma::bg.adjust.gcrma"
#  @seeclass
# }
#*/###########################################################################
setMethodS3("gcrmaSummary", "AffymetrixCelSet", function(this, path=NULL, name="gcrma", bgPath=NULL, normPath=NULL, summaryMethod="rlm", probePath=NULL, affinities=NULL, type="fullmodel",  indicesNegativeControl=NULL, opticalAdjust=TRUE, gsbAdjust=TRUE, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15*fast + (1-fast), fast=TRUE, ..., verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /gcrma/<data set name>/chip_data/<chip type>/
    path <- file.path(name, getName(this), "chip_data", getChipType(cdf));
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot compute expression summaries. Argument 'path' refers to the same path as the path of the raw probe level data: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # do the background correction step

  verbose && enter(verbose, "Performing background correction");

  dsBG <- bgAdjustGcrma(this, path=bgPath, name=name, probePath=probePath, affinities=affinities, type=type,  indicesNegativeControl=indicesNegativeControl, opticalAdjust=opticalAdjust, gsbAdjust=gsbAdjust, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15*fast + (1-fast), fast=fast, ..., verbose=verbose);

  verbose && exit(verbose);

  # normalisation step
  
  verbose && enter(verbose, "Normalising");
                  
  dsQN <- normalizeQuantile(dsBG, path=normPath, name=name, typesToUpdate="pm", ..., verbose=verbose);

  verbose && exit(verbose);

  verbose && enter(verbose, "Computing expression summaries");
  
  if (summaryMethod=="rlm") {
    rmaPlm <- RmaPlm(dsQN);
    fitResult <- fit(rmaPlm, verbose=verbose);
  } else {
    throw("Only rlm is currently supported as a summary method");
  }

  verbose && exit(verbose);
  
  res <- AffymetrixCelSet$fromFiles(path=getPath(rmaPlm), pattern="cel$");
  return(res);

  
})


###########################################################################/**
# @RdocMethod rmaSummary
#
# @title "Calculates the RMA expression summary"
#
# \description{
#  @get "title".
#
# Models the observed signal as normal background + exponential signal, then
# applies quantile normalisation followed by a summary by either median
# polish or robust regression.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path to save the expression summaries.}
#   \item{name}{Name of the data set containing the expression summaries.}
#   \item{bgPath}{Directory in which to store background-adjusted signals.}
#   \item{normPath}{Directory in which to store normalised signals.}
#   \item{summaryMethod}{Either "medianpolish" or "rlm".}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet" containing the expression summaries.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#
# \seealso{
#  See package \pkg{affy}.
#  @seeclass
# }
#*/###########################################################################
setMethodS3("rmaSummary", "AffymetrixCelSet", function(this, path=NULL, name="rma", bgPath=NULL, normPath=NULL, summaryMethod="rlm", ..., verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /rma/<data set name>/chip_data/<chip type>/
    path <- file.path(name, getName(this), "chip_data", getChipType(cdf));
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot compute expression summaries. Argument 'path' refers to the same path as the path of the raw probe level data: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # do the background correction step

  verbose && enter(verbose, "Performing background correction");

  dsBG <- bgAdjustRma(this, path=bgPath, name=name, ..., verbose=verbose);

  verbose && exit(verbose);

  verbose && enter(verbose, "Normalising");
  
  dsQN <- normalizeQuantile(dsBG, path=normPath, name=name, typesToUpdate="pm", ..., verbose=verbose);

  verbose && exit(verbose);

  verbose && enter(verbose, "Computing expression summaries");
  
  if (summaryMethod=="rlm") {
    rmaPlm <- RmaPlm(dsQN);
    fitResult <- fit(rmaPlm, verbose=verbose);
  } else {
    throw("Only rlm is currently supported as a summary method");
  }

  verbose && exit(verbose);
  
  res <- AffymetrixCelSet$fromFiles(path=getPath(rmaPlm), pattern="cel$");
  return(res);

  
})



setMethodS3("getFullName", "AffymetrixCelSet", function(this, parent=1, ...) {
  NextMethod("getFullName", this, parent=parent, ...);
})

setMethodS3("getName", "AffymetrixCelSet", function(this, ...) {
  NextMethod("getName", this, ...);
})





############################################################################
# HISTORY:
# 2006-12-01
# o Now as.character() reports the range of CEL header timestamps.
# 2006-11-07
# o Now getAverageFile() uses rowMedians() of R.native if available, 
#   otherwise a local version utilizing apply(). Same for rowMads().
# 2006-10-24
# o Added getAverageLog() and getAverageAsinh().
# o Added transforms and anti-transforms g() and h() to getAverageFile().
# o Changed the defaults from mean to median, and sd to mad for 
#   getAverageFile().
# o Added Rdoc comments to getAverageFile().
# 2006-10-10
# o Renamed rma and gcrma to rmaSummary and gcrmaSummary, to avoid clash
#   with existing functions.
# o Added gcrma() wrapper function.
# o Added rma() wrapper function.
# o Fixed a bug in getData() - default for argument "fields" contained "xy",
#   which is not a valid field (x, y are separate).
# 2006-10-02
# o Added getData().  Now getIntensities() works again and is just a wrapper
#   to getData().
# 2006-09-18
# o Now references to all requested average files are cached so it can
#   return the same object instead of creating a new one each time.
# 2006-09-16
# o Added getSiblings() to easily get other data sets for the same
#   samples.
# 2006-09-14
# o Added a read-buffer cache to readUnits() and getUnitIntensities().
# 2006-08-27
# o Added getAverageFile().
# 2006-08-26
# o Now getName() of a CEL set is inferred from the pathname:
#     path/to/<name>/chip_files/<"chip type">/
# 2006-08-21
# o Now AffymetrixCelSet inherits from AffymetrixFileSet.
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
