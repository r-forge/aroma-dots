###########################################################################/**
# @RdocClass AffymetrixCdfFile
#
# @title "The AffymetrixCdfFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCdfFile object represents a generic Affymetrix CDF file.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AffymetrixCdfFile", function(...) {
  this <- extend(AffymetrixFile(...), "AffymetrixCdfFile",
    "cached:.header" = NULL,
    "cached:.unitNames" = NULL,
    "cached:.unitSizes" = NULL,
    "cached:.cellIndices" = NULL,
    "cached:.isPm" = NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("clearCache", "AffymetrixCdfFile", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".header", ".unitNames", ".unitSizes", ".cellIndices", ".isPm")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("as.character", "AffymetrixCdfFile", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Filename: %s", getFilename(this)));
  s <- c(s, sprintf("Filesize: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s <- c(s, sprintf("Dimension: %s", paste(getDimension(this), collapse="x")));
  s <- c(s, sprintf("Number of cells: %d", nbrOfCells(this)));
  # Requires reading of data:
#  nbrOfPms <- sum(isPm(this));
#  s <- c(s, sprintf("Number of PM cells: %d (%.2f%%)", nbrOfPms, 100*nbrOfPms/nbrOfCells(this)));
  s <- c(s, sprintf("Number of units: %d", nbrOfUnits(this)));
  s <- c(s, sprintf("Cells per unit: %.2f", nbrOfCells(this)/nbrOfUnits(this)));
  # Requires that unit names are read:
#  s <- c(s, sprintf("Number of AFFX- units: %d", length(indexOf(this, "^AFFX-"))));
  s <- c(s, sprintf("Number of QC units: %d", nbrOfQcUnits(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixCdfFile object from a CDF file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of to the file.}
#  \item{path}{The path to the file.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an instance of @see "AffymetrixCdfFile" or its subclasses.
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# @author
#
# \seealso{
#   @seemethod "fromChipType".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "AffymetrixCdfFile", function(static, filename, path=NULL, ...) {
  # Arguments 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  # Assert that it is a CDF file
  header <- readCdfHeader(pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Try to define an instance of a subclass traversing bottom up.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  clazz <- Class$forName(class(static)[1]);
  for (className in rev(getKnownSubclasses(clazz))) {
    clazz <- Class$forName(className);
    tryCatch({
      res <- newInstance(clazz, pathname);
      return(res);
    }, error = function(ex) {})
  }

  newInstance(static, pathname);
}, static=TRUE)




###########################################################################/**
# @RdocMethod fromChipType
#
# @title "Defines an AffymetrixCdfFile object by chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.  
#  If a matching CDF file was not found, an error is thrown.
# }
#
# @author
#
# \seealso{
#   @seemethod "fromFile".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromChipType", "AffymetrixCdfFile", function(static, chipType, ...) {
  pathname <- static$findByChipType(chipType);
  if (is.null(pathname)) {
    throw("Could not create ", class(static)[1], " object. No CDF file with that chip type found: ", chipType);
  }

  fromFile(static, filename=pathname, path=NULL, ...);
}, static=TRUE)




###########################################################################/**
# @RdocMethod findByChipType
#
# @title "Locates a CDF file from its chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{...}{Not used.}
#  \item{.useAffxparser}{If @TRUE, @see "affxparser::findCdf" is used if
#    the CDF could not be located.}
# }
#
# \value{
#  Returns a pathname as a @character string to the first CDF file found.
#  If non CDF with requested chip type was found, @NULL is returned.
# }
#
# @author
#
# \details{
# }
#
# \seealso{
#   @seemethod "fromChipType".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("findByChipType", "AffymetrixCdfFile", function(static, chipType, ..., .useAffxparser=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Handle deprecated <chipType>-monocell CDFs specially
  pattern <- "-monocell$";
  if (regexpr(pattern, chipType) != -1) {
    newChipType <- gsub("-monocell$", ",monocell", chipType);
    parentChipType <- gsub("-monocell$", "", chipType);  # Remove tags

    # First, see if there is a new monocell, then use that
    pattern <- paste("^", newChipType, "[.](c|C)(d|D)(f|F)$", sep="");
    pathname <- findAnnotationDataByChipType(parentChipType, pattern);

    # Second, see if the old-named monocell is there
    if (is.null(pathname)) {
      pattern <- paste("^", chipType, "[.](c|C)(d|D)(f|F)$", sep="");
      pathname <- findAnnotationDataByChipType(parentChipType, pattern=pattern);
      if (!is.null(pathname)) {
        msg <- paste("Deprecated filename of monocell CDF detected. Rename CDF file by replacing dash ('-') with a comma (','): ", pathname, sep="");
        warning(msg);
      }
    }

    return(pathname);
  }

  pattern <- paste("^", chipType, "[.](c|C)(d|D)(f|F)$", sep="");
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- paste("^", chipType, "[.](c|C)(d|D)(f|F)[.]lnk$", sep="");
    args <- list(chipType=chipType, ...);
    args$pattern <- pattern;
    pathname <- do.call("findAnnotationDataByChipType", args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- filePath(pathname, expandLinks="any");
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the header of the CDF file"
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
#  Returns a @list structure as returned by @see "affxparser::readCdfHeader".
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
setMethodS3("getHeader", "AffymetrixCdfFile", function(this, ...) {
  if (is.null(header <- this$.header))
    header <- this$.header <- readCdfHeader(this$.pathname);
  header;
}, private=TRUE)


setMethodS3("getChipType", "AffymetrixCdfFile", function(this, fullname=TRUE, ...) {
  chipType <- getHeader(this)$chiptype;

  # Get the main chip type?
  if (!fullname) {
    # Handle '-monocell' specially
    pattern <- "^(.*)-(monocell)$";
    if (regexpr(pattern, chipType) != -1) {
      chipType <- gsub(pattern, "\\1", chipType);
      tags <- "monocell";
    } else {    
      pattern <- "^([^,]*)[,](.*)$";
      tags <- gsub(pattern, "\\2", chipType);
      tags <- strsplit(tags, split=",")[[1]];
      chipType <- gsub(pattern, "\\1", chipType);
    }
    attr(chipType, "tags") <- tags;
  }

  chipType;
})

setMethodS3("getDimension", "AffymetrixCdfFile", function(this, ...) {
  header <- getHeader(this);
  c(header$rows, header$cols);
})

setMethodS3("nbrOfRows", "AffymetrixCdfFile", function(this, ...) {
  as.integer(getDimension(this, ...)[1]);
})

setMethodS3("nbrOfColumns", "AffymetrixCdfFile", function(this, ...) {
  as.integer(getDimension(this, ...)[2]);
})

setMethodS3("nbrOfCells", "AffymetrixCdfFile", function(this, ...) {
  as.integer(prod(getDimension(this, ...)));
})

setMethodS3("nbrOfUnits", "AffymetrixCdfFile", function(this, ...) {
  getHeader(this)$probesets;
})

setMethodS3("nbrOfQcUnits", "AffymetrixCdfFile", function(this, ...) {
  getHeader(this)$qcprobesets;
})


###########################################################################/**
# @RdocMethod getUnitNames
#
# @title "Gets the names of each unit"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @vector of @character names.
# }
#
# \details{
#   Once read from file, this information is cached in memory for efficiency.
#   The cache can be cleared by calling \code{gc(cdf)}.
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
setMethodS3("getUnitNames", "AffymetrixCdfFile", function(this, units=NULL, ...) {
  if (is.null(names <- this$.unitNames))
    names <- this$.unitNames <- readCdfUnitNames(this$.pathname, ...);
  if (!is.null(units))
    names <- names[units];
  names;
})


###########################################################################/**
# @RdocMethod getUnitSizes
#
# @title "Gets the number of groups in each unit"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @vector of @integers.
# }
#
# \details{
#   Once read from file, this information is cached in memory for efficiency.
#   The cache can be cleared by calling \code{gc(cdf)}.
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
setMethodS3("getUnitSizes", "AffymetrixCdfFile", function(this, units=NULL, ...) {
  sizes <- this$.unitSizes;
  if (is.null(sizes)) {
    sizes <- readCdfGroupNames(this$.pathname);
    sizes <- restruct(this, sizes);
    sizes <- lapply(sizes, FUN=length);
    sizes <- unlist(sizes, use.names=FALSE);
    this$.unitSizes <- sizes;
  }

  if (!is.null(units))
    sizes <- sizes[units];

#  gc <- gc();

  sizes;
}, private=TRUE)


## setMethodS3("getUnitGroupSizes", "AffymetrixCdfFile", function(this, units=NULL, force=FALSE, ...) {
##   sizes <- this$.unitGroupSizes;
##   if (force || is.null(sizes)) {
##     sizes <- readCdfNbrOfCellsPerUnitGroup(this$.pathname);
##     sizes <- restruct(this, sizes);
##     this$.unitGroupSizes <- sizes;
##   }
## 
##   if (!is.null(units))
##     sizes <- sizes[units];
## 
##   sizes;
## }, private=TRUE)




###########################################################################/**
# @RdocMethod indexOf
#
# @title "Gets the indices of units by their names"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pattern}{A pattern to be used for identifying unit names of 
#      interest.  If @NULL, no regular expression matching is done.}
#   \item{names}{Names to be match exactly to the unit names.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @vector of @integers in [1,N] where N is the number of units
#  in this CDF structure.
# }
#
# @author
#
# \seealso{
#   @seemethod "getUnitNames".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("indexOf", "AffymetrixCdfFile", function(this, pattern=NULL, names=NULL, ...) {
  if (!is.null(names)) {
    idx <- match(names, getUnitNames(this));
  } else if (!is.null(pattern)) {
    idx <- grep(pattern, getUnitNames(this));
  } else {
    throw("Either argument 'names' or 'pattern' must be specified.");
  }

  idx;
})


###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Gets the cell indices unit by unit"
#
# \description{
#  @get "title" of all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Additional arguments passed to 
#      @see "affxparser::readCdfCellIndices".}
#   \item{force}{If @TRUE, cached values are ignored.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @list structure returned by 
#  @see "affxparser::readCdfCellIndices".
# }
#
# @author
#
# \seealso{
#   See @seemethod "setRestructor" to set a default re-constructor for
#   the returned CDF structure.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCellIndices", "AffymetrixCdfFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="getCellIndices", class=class(this)[1], 
                           chipType=getChipType(this), units=units, ...);
  id <- digest(key);
  res <- this$.cellIndices[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "getCellIndices.AffymetrixCdfFile(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read from CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading cell indices from CDF file");
  verbose && cat(verbose, "Pathname: ", this$.pathname);
  verbose && cat(verbose, "Units: ");
  verbose && str(verbose, units);
  verbose2 <- -as.integer(verbose)-1;
  cdf <- readCdfCellIndices(this$.pathname, units=units, ..., verbose=verbose);
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Restructuring");
  cdf <- restruct(this, cdf);  # Always call restruct() after a readCdfNnn()!
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (object.size(cdf) < 10e6) { # Cache only objects < 10MB.
    verbose && cat(verbose, "readUnits.AffymetrixCdfFile(): Updating cache");
    this$.cellIndices <- list();
    this$.cellIndices[[id]] <- cdf;
  }

  cdf;
})



setMethodS3("restruct", "AffymetrixCdfFile", function(this, cdf, ...) {
  # Rearrange CDF structure?
  fcn <- this$.restructor;
  if (!is.null(fcn))
    cdf <- fcn(cdf);
  cdf;
}, private=TRUE)





###########################################################################/**
# @RdocMethod setRestructor
#
# @title "Specifies a function through which"
#
# \description{
#  @get "title" of all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Additional arguments passed to 
#      @see "affxparser::readCdfCellIndices".}
# }
#
# \value{
#  Returns the @list structure returned by 
#  @see "affxparser::readCdfCellIndices".
# }
#
# \section{Requirements}{
#   The reconstructor function \emph{have to}:
#  \enumerate{
#   \item Accept a @list structure as its first argument.
#   \item Require no other arguments.
#   \item Return a @list structure of identical length as the input one.
#         In other words, it must not change the number of units nor
#         the order of units.
#  }
#
#  The reconstructor function \emph{may}:
#  \enumerate{
#   \item Rearrange the groups or change the number of groups, for
#         instance by utilizing @see "affxparser::applyCdfGroups".
#   \item Exclude some cells, for instance by excluding MMs.
#  }
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
setMethodS3("setRestructor", "AffymetrixCdfFile", function(this, fcn=NULL, ...) {
  if (is.null(fcn)) {
  } else if (is.function(fcn)) {
  } else {
    throw("Argument 'fcn' must be NULL or a function: ", mode(fcn));
  }
  if (!identical(this$.restructor, fcn)) {
    this$.restructor <- fcn;
    clearCache(this);
  }
  invisible(this);
}, private=TRUE)

setMethodS3("getRestructor", "AffymetrixCdfFile", function(this, ...) {
  this$.restructor;
}, private=TRUE)



###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads CDF data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probeset).
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Additional arguments passed to @see "affxparser::readCdfUnits".}
# }
#
# \value{
#  Returns the @list structure that @see "affxparser::readCdfUnits" returns
#  (possibly restructured).
# }
#
# \section{Caching}{
#   CDF data is neither cached in memory nor on file by this method.
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
# NOTE: getUnits() does not work because an S4 class stole it!!!
setMethodS3("readUnits", "AffymetrixCdfFile", function(this, units, ...) {
#  cdf <- readCdfUnits(this$.pathname, units=units, ...);
  cdf <- doCall("readCdfUnits", filename=this$.pathname, units=units, ...);
  restruct(this, cdf);  # Always call restruct() after a readCdfNnn()!
})




###########################################################################/**
# @RdocMethod isPm
#
# @title "Checks which cells (probes) are PMs and not"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Additional arguments passed to @see "affxparser::readCdfUnits".}
#   \item{force}{If @TRUE, cached results are ignored.}
#   \item{cache}{If @TRUE, results are cached.}
# }
#
# \value{
#  Returns a @logical @vector of length K, where K equals the total number
#  of cells in the requested units.  Note that the cells are ordered as
#  they occur in the units, that is, \emph{not} in incremental order.
# }
#
# \section{Caching}{
#   This method caches a @logical @vector of length N, when N equals the
#   number of cells on the array. The size of this vector is approximately
#   4*N bytes.  The vector indicates if a cell is a PM or not.
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
setMethodS3("isPm", "AffymetrixCdfFile", function(this, units=NULL, force=FALSE, cache=TRUE, ...) {
  isPm <- this$.isPm;
  if (force || is.null(isPm)) {
    if (cache) {
      # If caching, read all units
      cdf <- readCdfIsPm(this$.pathname);
      cdf <- restruct(this, cdf);  # Always call restruct() after a readCdfNnn()!
      isPm <- this$.isPm <- cdf;
    } else {
      # ...otherwise, read only a subset of units
      cdf <- readCdfIsPm(this$.pathname, units=units);
      cdf <- restruct(this, cdf);  # Always call restruct() after a readCdfNnn()!
      isPm <- cdf;
    }
  }

  if (cache && !is.null(units)) {
    isPm <- isPm[units];
  }

  # Return a vector
  isPm <- unlist(isPm, use.names=FALSE);

  isPm;
})


setMethodS3("identifyCells", "AffymetrixCdfFile", function(this, indices=NULL, from=1, to=nbrOfCells(this), types=c("all", "pmmm", "pm", "mm", "qc"), ..., sort=TRUE, .force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'types':
  if (is.null(types))
    types <- "all";

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Arguments 'from' and 'to':
  nbrOfCells <- nbrOfCells(this);
  from <- Arguments$getInteger(from, range=c(1, nbrOfCells));
  to <- Arguments$getInteger(to, range=c(1, nbrOfCells));

  # Argument 'indices':
  if (is.numeric(indices)) {
    getFraction <- (length(indices) == 1 && indices > 0 && indices < 1);
    if (getFraction) {
      by <- 1/indices;
    } else {
      indices <- Arguments$getIntegers(indices, range=c(1, nbrOfCells));
    }
  } else {
    getFraction <- FALSE;
  }

  if ("all" %in% types) {
    other <- 1:nbrOfCells;
  } else {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Check for cached results (already here)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Create a cache key (already here)
    verbose && enter(verbose, "Checking cache");
    chipType <- getChipType(this);
    key <- list(method="identifyCells", class=class(this)[1], chipType=chipType, indices=indices, from=from, to=to, types=types, sort=sort);
    comment <- sprintf("%s: %s", key$method, key$chipType);
    dirs <- c("aroma.affymetrix", chipType);
    if (!.force) {
      cache <- loadCache(key=key, dirs=dirs);
      if (!is.null(cache)) {
        verbose && exit(verbose, suffix="found cached value");
        return(cache);
      }
    }
    verbose && exit(verbose);
  }

  # Argument 'from':
  if (is.null(indices)) {
    indices <- seq(from=from, to=to, ...);
    indices <- as.integer(indices+0.5);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Intersect 'indices' and 'types'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!"all" %in% types) {
    verbose && enter(verbose, "Identifies cells of certain kind");
    indices <- unlist(getCellIndices(this), use.names=FALSE);
  
    other <- c();
    for (type in types) {
      if (type == "pm") {
        verbose && cat(verbose, "Using PM only");
        other <- c(other, indices[isPm(this)]);
      } else if (type == "mm") {
        verbose && cat(verbose, "Using MM only");
        other <- c(other, indices[!isPm(this)]);
      } else if (type == "pmmm") {
        verbose && cat(verbose, "Using PM & MM");
        other <- c(other, indices);
      } else if (type == "qc") {
        verbose && cat(verbose, "Using QC cells only");
        # Get cell indices for all non-regular units, i.e. QCs
        other <- c(other, setdiff(1:nbrOfCells, indices));
      }
    }

    other <- unique(other);
    verbose && exit(verbose);
  }

  if (is.null(indices)) {
    indices <- other;
    rm(other);
  } else {
    if (getFraction) {
      # Get the fraction from the already filtered cell indices
      indices <- other[seq(from=1, to=length(other), by=by)];
    } else if (!"all" %in% types) {
      indices <- intersect(indices, other);
    }
  }

  if (sort)
    indices <- sort(indices);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Save result to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!"all" %in% types) {
    saveCache(indices, key=key, comment=comment, dirs=dirs);
  }
  
  indices;
}, private=TRUE);


setMethodS3("getFirstCellIndices", "AffymetrixCdfFile", function(this, units=NULL, stratifyBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Trying to load cached results");
  chipType <- getChipType(this);
  key <- list(method="getFirstCellIndices", class=class(this)[1], chipType=chipType, stratifyBy=stratifyBy, restructor=body(this$.restructor));
  dirs <- c("aroma.affymetrix", chipType);
  res <- if (force) {
    NULL;
  } else {
    loadCache(key=key, dirs=dirs);
  }
  verbose && exit(verbose);

  if (is.null(res)) {
    verbose && enter(verbose, "Reading all cell indices (slow)");
    res <- getCellIndices(this, units=NULL, ..., stratifyBy=stratifyBy, verbose=verbose);
    verbose && exit(verbose);
      
    verbose && enter(verbose, "Extracting the first cell in each unit group");
    # For each unit and each group, get the index of the first cell.
    res <- applyCdfGroups(res, function(groups) {
      # For each group, pull out the first cell.
      lapply(groups, FUN=function(group) {
        # group$indices[1] == group[[1]][1] == ...
        list(indices=.subset(.subset2(group, 1), 1));
      })
    });
    verbose && exit(verbose);

    # Save to cache file
    verbose && enter(verbose, "Saving results to cache");
    saveCache(res, key=key, dirs=dirs);
    verbose && exit(verbose);
  }

  # Subset?
  if (!is.null(units))
    res <- res[units];

  res;
}, private=TRUE)


###########################################################################/**
# @RdocMethod compare
#
# @title "Checks if two AffymetrixCdfFile objects are equal"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{other}{The other @see "AffymetrixCdfFile" object to be compared to.}
#   \item{...}{Additional arguments passed to @see "affxparser::compareCdfs".}
# }
#
# \value{
#  Returns @TRUE if the two objects are equal, otherwise @FALSE.
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
setMethodS3("compare", "AffymetrixCdfFile", function(this, other, ...) {
  if (!inherits(other, "AffymetrixCdfFile"))
    return(FALSE);

  # Check if it is the same object
  if (equals(this, other))
    return(TRUE);

  res <- compareCdfs(getPathname(this), getPathname(other), ...);

  res;
}, private=TRUE)



###########################################################################/**
# @RdocMethod convert
#
# @title "Converts a CDF into the same CDF but with another format"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{chipType}{The chip type of the new CDF.}
#   \item{suffix}{A suffix added to the chip type of the new CDF.}
#   \item{sep}{A string separating the chip type and the suffix string.}
#   \item{path}{The path where to store the new CDF file.}
#   \item{...}{Additional arguments passed to @see "affxparser::convertCdf".}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the new CDF as an @see "AffymetrixCdfFile" object.
# }
#
# @author
#
# \seealso{
#   To compare two CDFs, see @seemethod "equals".
#   Internally @see "affxparser::convertCdf" is used.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("convert", "AffymetrixCdfFile", function(this, chipType=getChipType(this), suffix=NULL, sep="-", path="cdf", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Create the pathname of the destination CDF
  name <- paste(c(chipType, suffix), collapse=sep);
  dest <- sprintf("%s.cdf", name);
  dest <- Arguments$getWritablePathname(dest, path=path);

  # Convert CDF
  src <- getPathname(this);
  verbose2 <- -getThreshold(verbose);
  res <- convertCdf(src, dest, ..., verbose=verbose2);

  # Return an AffymetrixCdfFile object for the new CDF
  newInstance(this, dest);
})




###########################################################################/**
# @RdocMethod getGenomeInformation
#
# @title "Gets genome information for this chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{types}{A @character @vector specifying what type of genome 
#     information sets to search for.}
#   \item{...}{Not used.}
#   \item{force}{If @FALSE, cached information is retrieved, otherwise not.}
# }
#
# \value{
#  Returns a @see "GenomeInformation" object.
# }
#
# \examples{\dontrun{
#   @include "../incl/getGenomeInformation.Rex"
# }}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getGenomeInformation", "AffymetrixCdfFile", function(this, types=c("dChip"), ..., force=FALSE) {
  chipType <- getChipType(this, fullname=FALSE);

  gi <- this$.gi;
  if (is.null(gi) || force) {
    for (type in types) {
      tryCatch({
        if (type == "dChip") {
          gi <- DChipGenomeInformation$fromChipType(chipType);
        }
      }, error = function(ex) {})
    }
  
    if (is.null(gi)) {
      throw("Failed to retrieve genome information for this chip type: ", chipType);
    }

    this$.gi <- gi;
  }

  gi;
})


setMethodS3("getSnpInformation", "AffymetrixCdfFile", function(this, types=c("dChip"), ..., force=FALSE) {
  # Remove any suffices to get the "main" chip type.
  chipType <- getChipType(this, fullname=FALSE);

  si <- this$.si;
  if (is.null(si) || force) {
    for (type in types) {
      tryCatch({
        if (type == "dChip") {
          si <- DChipSnpInformation$fromChipType(chipType);
        }
      }, error = function(ex) {})
    }
  
    if (is.null(si)) {
      throw("Failed to retrieve SNP information for this chip type: ", chipType);
    }

    this$.si <- si;
  }

  si;
}, private=TRUE)


###########################################################################/**
# @RdocMethod convertUnits
#
# @title "Gets and validates unit indices"
#
# \description{
#  @get "title" either by unit names or by a unit indices (validation).
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{Either a @character @vector with unit names, or an @integer
#     @vector with unit indices to be validated.  
#     If @NULL, all unit indices are returned.}
#   \item{keepNULL}{If @TRUE, @NULL returns @NULL.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector with unit indices.
#  If some units are non existing, an error is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("convertUnits", "AffymetrixCdfFile", function(this, units=NULL, keepNULL=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
    # Return all units
    if (keepNULL)
      return(NULL);
    units <- 1:nbrOfUnits(this);
  } else if (is.character(units)) {
    # Identify units by their names
    unitNames <- units;
    units <- indexOf(this, names=unitNames);
    if (any(is.na(units))) {
      missing <- unitNames[is.na(units)];
      nmissing <- length(missing);
      if (nmissing > 10)
        missing <- c(missing[1:10], "...");
      throw("Argument 'units' contains ", nmissing, " unknown unit names: ", 
                                             paste(missing, collapse=", "));
    }
  } else {
    # Validate unit indices
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(this)));
  }

  units;
}, private=TRUE)


############################################################################
# HISTORY:
# 2007-03-26
# o Added a few more gc().
# o BUG FIX: isPm() did not work when querying a subset of the units.
# 2007-02-22
# o Now findByChipType() recognizes Windows shortcuts too.
# 2007-02-21
# o Now findByChipType() passes '...' to underlying function.
# 2007-02-14
# o BUG FIX: When "tagifying" monocell, getSnpInformation() and
#   getGenomeInformation() was looking for the incorrect chip type.
# 2007-02-12
# o Added argument 'main' to getChipType().
# 2007-02-08
# o Now findByChipType() handles monocell CDFs specially; monocell CDFs can
#   still be put in the same directory as the parent CDF.
# 2007-02-06
# o Added findByChipType().
# 2007-01-16
# o Now all cache keys contains method name, class name, and chip type.
# 2007-01-10
# o Reordered internally in createMonoCell() preparing for code to read 
#   *and* write monocell CDFs in chunks.  It should not be too hard.
#   We need to update affxparser with writeCdfHeader(), writeCdfQcUnits()
#   and writeCdfUnits(), and are basically already in there, but as
#   private functions.
# o Removed some unnecessary group fields in 'destUnits' saving us approx
#   220-230Mb out of 1.9GB for the Mapping250K_Nsp.  Still need to find
#   another solution to get down below 1GB. One thing that takes up a lot
#   of memory is that the unit and group directions are stored as character
#   strings and not integers.
# o BUG FIX: The most recent createMonoCell() would create CDFs with all
#   cell indices being equal to one.  Added more verbose output and some
#   garbage collection to this function too.
# 2007-01-06
# o Added argument 'force' to getCellIndices().
# o Now getCellIndices() only caches object < 10 MB RAM.
# o Optimized identifyCells() to only cache data in rare cases.
# 2006-12-18 /KS
# o Made global replacement "block" -> "group".
# 2006-12-14
# o Added convertUnits().
# 2006-09-27
# o Now fromFile() tries to create an instance of the subclasses (bottom up)
#   first.  This will make it possible to automatically define SNP CDFs.
# 2006-09-26
# o Now getGenomeInformation() and getSnpInformation() ignores suffices of
#   the chip-type string. This makes it possible to retrive annotation data
#   also for custom chips.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-16
# o Added getGenomeInformation() and stextChipType().
# 2006-09-14
# o BUG FIX: Fractional value of 'indices' of identifyCells().
# 2006-09-10
# o BUG FIX: createMonoCell() where resetting the cell indices for each
#   chunk.
# o Simple benchmarking of createMonoCell(): IBM Thinkpad A31 1.8GHz 1GB:
#   Mapping50K_Hind to mono cells CDF takes ~13 mins.  Again, it is 
#   writeCdf() that is slow.  KH is working on improving this.
# 2006-09-08
# o Added equals() to compare to CDF object.
# o Added convert() to convert a CDF into another version by convertCdf().
# 2006-08-25
# o Added getFirstCellIndices().  This may be used by methods to identify
#   unit or unit groups for which no probe signals have been assigned yet.
# 2006-08-24
# o Added reconstruct().
# o Added Rdoc comments.
# 2006-08-19
# o Added getUnitSizes().  Useful in order to store parameter estimates
#   of probeset-summary models.
# 2006-08-11
# o Created.
############################################################################
