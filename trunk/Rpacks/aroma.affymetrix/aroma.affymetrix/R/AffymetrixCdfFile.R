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


############################################################################
# HISTORY:
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
