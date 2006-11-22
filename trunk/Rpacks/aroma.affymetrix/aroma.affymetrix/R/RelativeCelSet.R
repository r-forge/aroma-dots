###########################################################################/**
# @RdocClass RelativeCelSet
#
# @title "The RelativeCelSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents an AffymetrixCelSet whose signals are calculated
#  relative to the signals of a reference AffymetrixCelFile.
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "AffymetrixCelSet".}
#  \item{...}{Arguments passed to the constructor of @see "AffymetrixCelSet".}
#  \item{reference}{A reference @see "AffymetrixCelFile" object with a
#    compatible chip type.}
#  \item{relativeFields}{A @character @vector specifying for which fields 
#    relative signals should be calculated.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("RelativeCelSet", function(dataSet=NULL, reference=NULL, relativeFields=c("intensities"), ...) {
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet")) {
      throw("Argument 'dataSet' is not a AffymetrixCelSet object: ", 
                                                           class(dataSet)[1]);
    }

    if (is.null(reference)) {
      throw("Argument 'reference' must not be NULL.");
    }
  }

  if (!is.null(reference)) {
    if (!inherits(reference, "AffymetrixCelFile")) {
      throw("Argument 'reference' is not a AffymetrixCelFile object: ", 
                                                         class(reference)[1]);
    }
  }

  if (!is.null(dataSet)) {
    cdf <- getCdf(dataSet);
    refCdf <- getCdf(reference);
    if (!equals(cdf, refCdf)) {
      throw("The CDFs for the data set and the reference file are not compatible: ", getChipType(cdf), " != ", getChipTYpe(refCdf));
    }
  }

  extend(AffymetrixCelSet(files=as.list(dataSet), ...), "RelativeCelSet",
    .reference = reference,
    .relativeFields = relativeFields 
  )
})

setMethodS3("getReferenceFile", "RelativeCelSet", function(this, ...) {
  this$.reference;
})

setMethodS3("getRelativeFields", "RelativeCelSet", function(this, ...) {
  this$.relativeFields;
})

setMethodS3("getData", "RelativeCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data from the data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- NextMethod("getData", this, ..., verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate ratios for relative fields
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  relFields <- intersect(names(data), getRelativeFields(this));
  if (length(relFields) > 0) {
    # Retrieve data from the reference data file
    args <- list(getReferenceFile(this), ...);
    args <- args[names(args) != "fields"];
    args <- c(args, list(fields=relFields, verbose=verbose));
    relData <- do.call("getData", args);

    for (field in relFields)
      data[[field]] <- data[[field]] / relData[[field]];
  }

  data;
})


setMethodS3("readUnits", "RelativeCelSet", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
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
  if (!is.list(units)) {
    # Always ask for CDF information from the CDF object!
    verbose && enter(verbose, "Retrieving CDF unit information");
    units <- readUnits(getCdf(this), units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }
  data <- readCelUnits(pathnames, cdf=units, ...);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate the relative signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  firstGroup <- data[[1]][[1]];
  relFields <- intersect(names(firstGroup), getRelativeFields(this));
  if (length(relFields) > 0) {
    # Retrieve data from the reference data file
    pathname <- getPathname(getRelativeFile(this));
    relData <- readCelUnits(pathname, cdf=units, ...);

    # For each unit
    for (uu in seq(along=data)) {
      dataUU <- .subset2(data, uu);
      relDataUU <- .subset2(relData, uu);
      # For each group
      for (gg in seq(along=dataUU)) {
        dataGG <- .subset2(dataUU, gg);
        relDataGG <- .subset2(relDataUU, gg);
        # For each "relative" fields
        for (field in relFields) {
          dataGG[[field]] <- dataGG[[field]] / relDataGG[[field]];
        }
        dataUU[[gg]] <- dataGG;
      }
      data[[uu]] <- dataUU;
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  this$.readUnitsCache <- list()
  this$.readUnitsCache[[key]] <- data
  verbose && cat(verbose, "readUnits(): Updated cache");

  verbose && exit(verbose);

  data;
})




############################################################################
# HISTORY:
# 2006-11-02
# o Created.
############################################################################
