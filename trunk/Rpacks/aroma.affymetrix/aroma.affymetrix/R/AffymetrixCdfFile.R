				
#
# @title "The AffymetrixCdfFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCdfFile object represents a generic Affymetrix 
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
  extend(AffymetrixFile(...), "AffymetrixCdfFile",
    "cached:.header" = NULL,
    "cached:.unitNames" = NULL,
    "cached:.unitSizes" = NULL,
    "cached:.isPm" = NULL
  )
})


setMethodS3("as.character", "AffymetrixCdfFile", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Filename: %s", getFilename(this)));
  s <- c(s, sprintf("Filesize: %.2fMb", getFileSize(this)/1024^2));
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
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})



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
#  Returns an @see "AffymetrixCdfFile" object.  
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

  newInstance(static, pathname);
}, static=TRUE)




###########################################################################/**
# @RdocMethod fromChipType
#
# @title "Defines an AffymetrixCdfFile object by chip type"
#
# \description{
#  @get "title" utilizing the @see "affxparser::findCdf" method.
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
  pathname <- findCdf(chipType);
  if (is.null(pathname)) {
    throw("Could not create ", class(static), " object. No CDF file with that chip type found: ", chipType);
  }

  fromFile(static, filename=pathname, path=NULL, ...);
}, static=TRUE)




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
})

setMethodS3("getChipType", "AffymetrixCdfFile", function(this, ...) {
  getHeader(this)$chiptype;
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

  sizes;
})



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
setMethodS3("getCellIndices", "AffymetrixCdfFile", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Reading cell indices from CDF file");
  cdf <- readCdfCellIndices(this$.pathname, units=units, ...);
  verbose && exit(verbose);

  verbose && enter(verbose, "Restructuring");
  cdf <- restruct(this, cdf);  # Always call restruct() after a readCdfNnn()!
  verbose && exit(verbose);
  cdf;
})






setMethodS3("restruct", "AffymetrixCdfFile", function(this, cdf, ...) {
  # Rearrange CDF structure?
  fcn <- this$.restructor;
  if (!is.null(fcn))
    cdf <- fcn(cdf);
  cdf;
}, protected=TRUE)





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
})

setMethodS3("getRestructor", "AffymetrixCdfFile", function(this, ...) {
  this$.restructor;
}, protected=TRUE)


# NOTE: getUnits() does not work because an S4 class stole it!!!
setMethodS3("readUnits", "AffymetrixCdfFile", function(this, units, ...) {
  cdf <- readCdfUnits(this$.pathname, units=units, ...);
  restruct(this, cdf);  # Always call restruct() after a readCdfNnn()!
})



setMethodS3("isPm", "AffymetrixCdfFile", function(this, units=NULL, ...) {
  isPm <- this$.isPm;
  if (is.null(isPm)) {
    cdf <- readCdfIsPm(this$.pathname);
    cdf <- restruct(this, cdf);  # Always call restruct() after a readCdfNnn()!
    isPm <- this$.isPm <- cdf;
  }

  if (!is.null(units)) {
    isPm <- isPm[units];
  }

  isPm <- unlist(isPm, use.names=FALSE);
  isPm;
})


setMethodS3("identifyCells", "AffymetrixCdfFile", function(this, indices=NULL, from=1, to=nbrOfCells(this), ..., types=c("all", "pmmm", "pm", "mm", "qc"), sort=TRUE, .force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'type':
  if (is.null(types))
    types <- "all";

  nbrOfCells <- nbrOfCells(this);
  from <- Arguments$getInteger(from, range=c(1, nbrOfCells));
  to <- Arguments$getInteger(to, range=c(1, nbrOfCells));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Create a cache key (already here)
  key <- list(method=sprintf("identifyCells(<%s>)", class(this)[1]), chipType=getChipType(this), indices=indices, from=from, to=to, types=types, sort=sort);
  comment <- sprintf("%s: %s", key$method, key$chipType);

  if (!.force) {
    cache <- loadCache(key=key);
    if (!is.null(cache))
      return(cache);
  }

  # Argument 'indices':
  if (is.numeric(indices)) {
    getFraction <- (length(indices) == 1 && indices > 0 && indices < 1);
    if (!getFraction) {
      indices <- Arguments$getIntegers(indices, range=c(1, nbrOfCells));
    }
  } else {
    getFraction <- FALSE;
  }

  # Argument 'from':
  if (is.null(indices)) {
    indices <- seq(from=from, to=to, ...);
    indices <- as.integer(indices+0.5);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  if ("all" %in% types)
    return(indices);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Intersect 'indices' and 'type'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
  if (is.null(indices)) {
    indices <- other;
  } else {
    if (getFraction) {
      # Get the fraction from the already filtered cell indices
      indices <- other[seq(from=1, to=length(other), by=1/indices)];
    } else {
      indices <- intersect(indices, other);
    }
  }

  if (sort)
    indices <- sort(indices);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Save result to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  saveCache(indices, key=key, comment=comment);
  
  indices;
}, protected=TRUE);


setMethodS3("getFirstCellIndices", "AffymetrixCdfFile", function(this, units=NULL, stratifyBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Trying to load cached results");
  key <- list(method="getFirstCellIndices.AffymetrixCdfFile", chipType=getChipType(this), stratifyBy=stratifyBy, restructor=body(this$.restructor));
  res <- if (force) NULL else loadCache(key=key);
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
    saveCache(key=key, res);
    verbose && exit(verbose);
  }

  # Subset?
  if (!is.null(units))
    res <- res[units];

  res;
}, protected=TRUE)


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
#  Returns @TRUE if the two objects are equal, otherwise @FALSE.}
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
})

# To be removed.
setMethodS3("equals", "AffymetrixCdfFile", function(...) {
  NextMethod("equals", ...);
})


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
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

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





setMethodS3("createMonoCell", "AffymetrixCdfFile", function(this, chipType=getChipType(this), suffix="monocell", sep="-", path="cdf", ..., nbrOfCellsPerField=1, nbrOfUnitsPerChunks=5000, verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rearrangeCells <- function(units, offset=0, hasBlocks=TRUE, ...) {
    rearrangeBlock <- function(block, idxs, ...) {
      y = (idxs-1) %/% ncols;
      x = (idxs-1) - ncols*y;
      block$y <- y;
      block$x <- x;
      block$indices <- idxs;
      block;
    } # rearrangeBlock()

    nbrOfCells <- lapply(units, FUN=function(unit) unit$ncells);
    nbrOfCells <- sum(unlist(nbrOfCells, use.names=FALSE));

    cells <- seq(from=offset+1, to=offset+nbrOfCells);

    verbose && printf(verbose, "Units: ");
    if (hasBlocks) {
      for (kk in seq(along=units)) {
        if (verbose) {
          if (kk %% 1000 == 0) {
            printf(verbose, "%d, ", kk);
          } else if (kk %% 100 == 0) {
            cat(".");
          }
        }
        # blocks <- units[[kk]]$blocks;
        blocks <- .subset2(.subset2(units, kk), "blocks");
        for (ll in seq(along=blocks)) {
          block <- .subset2(blocks, ll);
          # Number of cells in this block
          # nindices <- length(block$indices);
          nindices <- length(.subset2(block, "indices"));
          head <- 1:nindices;
          # idxs <- cells[head];
          idxs <- .subset(cells, head);
          # cells <- cells[<tail>];
          cells <- .subset(cells, (nindices+1):length(cells));
          blocks[[ll]] <- rearrangeBlock(block, idxs);
        }
        units[[kk]]$blocks <- blocks;
      }
    } else {
      for (kk in seq(along=units)) {
        if (verbose) {
          if (kk %% 1000 == 0) {
            printf(verbose, "%d, ", kk);
          } else if (kk %% 100 == 0) {
            cat(".");
          }
        }
        # block <- units[[kk]]
        block <- .subset2(units, kk);
        # Number of cells in this block
        # nindices <- length(block$indices);
        nindices <- length(.subset2(block, "indices"));
        head <- 1:nindices;
        # idxs <- cells[head];
        idxs <- .subset(cells, head);
        # cells <- cells[<tail>];
        cells <- .subset(cells, (nindices+1):length(cells));
        block <- rearrangeBlock(block, idxs);
        units[[kk]] <- block;
      }
    }
    verbose && printf(verbose, "\n");

    units;
  } # rearrangeCells()
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'nbrOfCellsPerField':
  nbrOfCellsPerField <- Arguments$getIndices(nbrOfCellsPerField);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the pathname of the source
  src <- getPathname(this);
  src <- Arguments$getReadablePathname(src);

  # Create the pathname of the destination CDF
  name <- paste(c(chipType, suffix), collapse=sep);
  dest <- sprintf("%s.cdf", name);
  dest <- Arguments$getWritablePathname(dest, path=path, mustNotExist=TRUE);

  # Assure source and destination is not identical
  if (identical(src, dest)) {
    throw("Cannot not create CDF file. Destination is same as source: ", src);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fields to be kept
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Number of cells to keep in each block field
  fidx <- 1:nbrOfCellsPerField;

  verbose && printf(verbose, "Number of cells per block field: %d\n", 
                                                       nbrOfCellsPerField);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # QC units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Keep all QC units  
  destQcUnits <- readCdfQc(src);
  nbrOfQcUnits <- length(destQcUnits);
  nbrOfQcCells <- lapply(destQcUnits, FUN=function(unit) unit$ncells);
  nbrOfQcCells <- sum(unlist(nbrOfQcCells, use.names=FALSE));
  verbose && printf(verbose, "Number of QC cells: %d in %d QC units\n", 
                                                nbrOfQcCells, nbrOfQcUnits);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chip layout
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unitSizes <- readCdfNbrOfCellsPerUnitGroup(src);
  unitSizes <- unlist(lapply(unitSizes, FUN=length), use.names=FALSE);
  unitSizes <- nbrOfCellsPerField * unitSizes;
  nbrOfCells <- sum(unitSizes);

  totalNbrOfCells <- nbrOfCells + nbrOfQcCells;
  verbose && printf(verbose, "Total number of cells: %d\n", totalNbrOfCells);

  # Figure out a best fit square layout of this number of cells
  side <- as.integer(floor(sqrt(totalNbrOfCells)));
  nrows <- ncols <- side;
  if (nrows*ncols < totalNbrOfCells) {
    nrows <- as.integer(nrows + 1);
    if (nrows*ncols < totalNbrOfCells) {
      ncols <- as.integer(ncols + 1);
    }
  }
  verbose && printf(verbose, "Best array dimension: %dx%d (=%d cells, i.e. %d left-over cells)\n", nrows, ncols, nrows*ncols, nrows*ncols - totalNbrOfCells);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create all new units in chunks
  nbrOfUnits <- nbrOfUnits(this);

  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

  unitsToDo <- 1:nbrOfUnits;

  destUnits <- vector("list", nbrOfUnits);
  nbrOfChunks <- ceiling(nbrOfUnits / nbrOfUnitsPerChunks);

  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n", 
                                         nbrOfChunks, nbrOfUnitsPerChunks);
  head <- 1:nbrOfUnitsPerChunks;
  verbose && enter(verbose, "Extracting unit data");
  fields <- c("x", "y", "indices", "pbase", "tbase", "atom", "indexpos");
  fields <- c("pbase", "tbase", "atom", "indexpos");
  count <- 1;
  idxOffset <- as.integer(0);
  while (length(unitsToDo) > 0) {
    if (length(unitsToDo) < nbrOfUnitsPerChunks) {
      head <- 1:length(unitsToDo);
    }
    units <- unitsToDo[head];
    verbose && printf(verbose, "Chunk #%d of %d (%d units)\n", 
                                        count, nbrOfChunks, length(units));

    unitsToDo <- unitsToDo[-head];

    # Read CDF structure
    srcUnits <- readCdf(src, units=units);

    srcUnits <- lapply(srcUnits, function(unit) {
      unit$blocks <- lapply(unit$blocks, function(block) {
        block[fields] <- lapply(.subset(block, fields), FUN=.subset, fidx);
        block$natoms <- nbrOfCellsPerField;
        block$ncellsperatom <- 1;
        idxs <- idxOffset + 1:nbrOfCellsPerField;
        y <- (idxs-1) %/% ncols;
        block$y <- y;
        block$x <- (idxs-1) - ncols*y;
        block$indices <- idxs;
        idxOffset <<- idxOffset + nbrOfCellsPerField;
        block;
      })
      ncells <- length(unit$blocks)*nbrOfCellsPerField;
      unit$ncells <- ncells;
      unit$natoms <- ncells;
      unit$ncellsperatom <- 1;
      unit;
    })

    # Store the transformed units
    destUnits[units] <- srcUnits;
    names(destUnits)[units] <- names(srcUnits);

    rm(srcUnits, units); # Not needed anymore

    count <- count + 1;
  } # while (length(unitsToDo) > 0)
  rm(unitsToDo, head, fields, fidx, count);
  nbrOfCells <- lapply(destUnits, FUN=function(unit) unit$ncells);
  nbrOfCells <- sum(unlist(nbrOfCells, use.names=FALSE));
  verbose && printf(verbose, "Number unit cells extracted: %d\n", nbrOfCells);
  verbose && exit(verbose);

  verbose && enter(verbose, "Rearranging QC unit cell indices");
  destQcUnits <- rearrangeCells(destQcUnits, offset=nbrOfCells, hasBlocks=FALSE, verbose=verbose);
  verbose && enter(verbose, "Validating QC unit cell indices");
  cells <- unlist(lapply(destQcUnits, FUN=function(unit) unit$indices), use.names=FALSE);
  udcells <- unique(diff(cells));
  if (!identical(udcells, 1:1)) {
    throw("Internal CDF error: The cell indices for the QC units are not contiguous: ", paste(udcells, collapse=", "));
  }
  verbose && exit(verbose);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # CDF header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creaging CDF header with source CDF as template");
  destHeader <- readCdfHeader(src);
  destHeader$nrows <- nrows;
  destHeader$ncols <- ncols;
  verbose && exit(verbose);

  # Write new CDF
  verbose && enter(verbose, "Writing new CDF to file");
  verbose && printf(verbose, "Pathname: %s\n", dest);
  verbose2 <- as.integer(verbose)-1;
  verbose2 <- 2;
  writeCdf(dest, cdfheader=destHeader, cdf=destUnits, 
                                   cdfqc=destQcUnits, ..., verbose=verbose2);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Verifying the CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Verifying the written CDF");
  # Checking header
  header <- readCdfHeader(dest);
  if ((header$nrows != nrows) || (header$ncols != ncols)) {
    throw(sprintf("Failed to create a valid mono-cell CDF: The dimension of the written CDF does not match the intended one: (%d,%d) != (%d,%d)", header$nrows, header$ncols, nrows, ncols));
  }

  # Checking cell indices
  cells <- readCdfCellIndices(dest);
  cells <- unlist(cells, use.names=FALSE);
  if (!identical(unique(diff(cells)), 1:1)) {
    throw("Failed to create a valid mono-cell CDF: The cell indices are not contiguous: ", paste(udcells, collapse=", "));
  }
  verbose && exit(verbose);

  # Return an AffymetrixCdfFile object for the new CDF
  newInstance(this, dest);
})



############################################################################
# HISTORY:
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
