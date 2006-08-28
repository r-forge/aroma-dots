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
# @examples "../incl/AffymetrixCelSet.Rex"
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
    "cached:.fileSize" = NULL
  )
})


setMethodS3("clone", "AffymetrixCelSet", function(this, ..., verbose=TRUE) {
  # Clone itself and the files.  The call below will clear the cache!
  object <- NextMethod("clone", clear=FALSE, ...);

  # Clone the CDF (this will update the CDF of all file object)
  cdf <- clone(getCdf(object));
  setCdf(object, cdf);

  object;
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
setMethodS3("setCdf", "AffymetrixCelSet", function(this, cdf, ...) {
  # Nothing to do?
  oldCdf <- getCdf(this);
  if (equals(cdf, oldCdf))
    return(invisible(this));

  # Set the CDF for all CEL files
  lapply(this, setCdf, cdf, ...);

  # Have to clear the cache 
  clearCache(this);

  invisible(this);
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
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("Number of arrays: %d", nbrOfArrays(this)));
  s <- c(s, sprintf("Total file size: %.2fMb", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("fromFiles", "AffymetrixCelSet", function(static, path="cel/", pattern="[.](c|C)(e|E)(l|L)$", ..., fileClass="AffymetrixCelFile") {
  this <- fromFiles.AffymetrixFileSet(static, path=path, pattern=pattern, ..., fileClass=fileClass);
  # Use the same CDF object for all CEL files.
  setCdf(this, getCdf(this));
  this;
})


setMethodS3("getName", "AffymetrixCelSet", function(this, ...) {
  # The name of a file set is inferred from the pathname of the directory
  # of the set, i.e. path/to/<name>/chip_files/<"chip type">/
  # Get the path of this file set
  path <- getPath(this);

  while(TRUE) {
    path <- dirname(path);
    # If that directory is named 'chip_data' (or 'chipdata') go up on level
    name <- basename(path);
    pattern <- "^chip[_]*data$";
    if (regexpr(pattern, name) == -1)
      break;
  }
  
  name;
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
setMethodS3("getIntensities", "AffymetrixCelSet", function(this, indices=NULL, verbose=FALSE, ...) {
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
  verbose && enter(verbose, "Getting cell signals for ", nbrOfArrays, " arrays.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading cell signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating the return matrix
  res <- matrix(NA, nrow=nbrOfCells, ncol=nbrOfArrays);

  for (kk in seq(length=nbrOfArrays)) {
    verbose && enter(verbose, "Array #", kk, " of ", nbrOfArrays);
    dataFile <- this$files[[kk]];
    res[,kk] <- getIntensities(dataFile, indices=indices, verbose=verbose);

    verbose && exit(verbose);
  }


  verbose && exit(verbose);

  res;
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

  res;
})



setMethodS3("readUnits", "AffymetrixCelSet", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the pathnames of all CEL files
  pathnames <- unlist(lapply(this, getPathname), use.names=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Cached values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.list(units)) {
    res <- readCelUnits(pathnames, cdf=units, ...);
  } else {
    # Always ask for CDF information from the CDF object!
    cdf <- readUnits(getCdf(this), units=units);
    res <- readCelUnits(pathnames, cdf=cdf, ...);
  }

  res;
})


setMethodS3("getAverageFile", "AffymetrixCelSet", function(this, name=NULL, indices="remaining", mean=c("mean", "median"), sd=c("sd", "mad"), na.rm=FALSE, ..., cellsPerChunk=moreCells*10^7/length(this), moreCells=1, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mean':
  if (is.character(mean)) {
    mean <- match.arg(mean);
    meanName <- mean;
    if (mean == "mean")
      mean <- base::rowMeans;
  } else if (is.function(mean)) {
    meanName <- "customMean";
  } else {
    throw("Argument 'mean' must be either a character or a function: ", mode(mean));
  }

  # Argument 'sd':
  if (is.character(sd)) {
    sd <- match.arg(sd);
    sdName <- sd;
    if (sd == "sd")
      sd <- rowSds;
  } else if (is.function(sd)) {
    sdName <- "customSd";
  } else {
    throw("Argument 'sd' must be either a character or a function: ", mode(sd));
  }

  # Argument 'name':
  if (is.null(name))
    name <- sprintf("average-%s-%s", meanName, sdName);

  # Argument 'indices':
  df <- as.list(this)[[1]];
  nbrOfCells <- getHeader(df)$total;
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create CEL file to store the average array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a private filename (with a dot prefix) to make sure it is not
  # identified as a regular CEL file when the directory is scanned for files.
  filename <- sprintf(".%s.CEL", name);
  res <- createFrom(df, filename=filename, path=getPath(this), verbose=verbose);
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

    verbose && enter(verbose, "Estimating averages and standard deviations");
    if (na.rm)
      n <- apply(X, MARGIN=1, FUN=function(x) { sum(!is.na(x)) });

    # Calculate the mean signal    
    mu <- mean(X, na.rm=na.rm);          # Special mean()!

    # Calculate the standard deviation of the signals
    sigma <- sd(X, mean=mu, na.rm=na.rm);   # Special sd()!
    verbose && exit(verbose);

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

  res;  
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

############################################################################
# HISTORY:
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
