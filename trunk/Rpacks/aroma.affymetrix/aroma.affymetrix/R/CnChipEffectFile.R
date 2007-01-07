###########################################################################/**
# @RdocClass CnChipEffectFile
#
# @title "The CnChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in a copy-number probe-level
#  models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpChipEffectFile".}
#   \item{combineAlleles}{A @logical indicating if the signals from allele A and
#      allele B are combined or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically part of a @see "CnChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("CnChipEffectFile", function(..., combineAlleles=FALSE) {
  extend(SnpChipEffectFile(...), "CnChipEffectFile",
    combineAlleles = combineAlleles
  )
})

setMethodS3("getParameters", "CnChipEffectFile", function(this, ...) {
  params <- NextMethod(generic="getParameters", object=this, ...);
  params$combineAlleles <- this$combineAlleles;
  params;
})


setMethodS3("getCellIndices", "CnChipEffectFile", function(this, ..., force=FALSE, .cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "getCellIndices.CnChipEffectFile()");



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force || .cache) {
    chipType <- getChipType(getCdf(this));
    params <- getParameters(this);
    key <- list(chipType=chipType, params=params, ...);
    key <- digest(key);
  }

  if (!force) {
    res <- this$.cellIndices[[key]];
    if (identical(res, "on-file")) {
      verbose && cat(verbose, "Cached on file");
      res <- loadCache(list(key));
    }
    if (!is.null(res)) {
      verbose && cat(verbose, "Returning cached value");
      verbose && exit(verbose);
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get and restructure cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cells <- NextMethod("getCellIndices", this, ..., force=force, 
                                          .cache=FALSE, verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Combine alleles?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If combining alleles, return only every second group.
  # In order to improve readability we merge the names of alleles groups
  # combined, e.g. groups 'C' and 'G' become group 'CG'.
  if (this$combineAlleles) {
    verbose && enter(verbose, "Combining alleles");
    # Hard-wiring 1, 2 & 4 groups speed up things 3 times!
    cells <- applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
      names <- names(groups);
      if (ngroups == 4) {
        groups <- .subset(groups, c(1,3));
        names <- paste(.subset(names, c(1,3)), .subset(names, c(2,4)), sep="");
      } else if (ngroups == 2) {
        groups <- .subset(groups, 1);
        names <- paste(.subset(names, 1), .subset(names, 2), sep="");
      } else if (ngroups == 1) {
        groups <- .subset(groups, 1);
      } else {
        groups <- .subset(groups, odds);
        odds <- seq(from=1, to=ngroups, by=2);
        evens <- seq(from=2, to=ngroups, by=2);
        names <- paste(.subset(names, odds), .subset(names, evens), sep="");
      }
      names(groups) <- names;
      groups;
    })
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.cache) {
    # In-memory or on-file cache?
    size <- object.size(cells);
      verbose && printf(verbose, "Object size: %.1fMB\n", size/1024^2);
    if (size < 10e6) { 
      # In-memory cache for objects < 10Mb.
      this$.cellIndices <- list();
      this$.cellIndices[[key]] <- cells;
      verbose && cat(verbose, "Result cached in memory");
    } else {
      # On-file cache
      # Keep, in-memory cache.
      if (!is.list(this$.cellIndices))
        this$.cellIndices <- list();
      saveCache(cells, key=list(key));
      this$.cellIndices[[key]] <- "on-file";
      verbose && cat(verbose, "Result cached to file");
    }
  }

  verbose && exit(verbose);

  cells;
})


setMethodS3("readUnits", "CnChipEffectFile", function(this, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Check for cached data
  key <- digest(list(class=class(this), combineAlleles=this$combineAlleles, ...));
  res <- this$.readUnitsCache[[key]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.CnChipEffectFile(): Returning cached data");
    return(res);
  }

  # Retrieve the data
  res <- NextMethod("readUnits", this, ..., force=TRUE, cache=FALSE, verbose=verbose);


  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.CnChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[key]] <- res;
  }

  res;
})


############################################################################
# HISTORY:
# 2007-01-06
# o Now getCellIndices() caches large objects to file and small in memory.
# o Made getCellIndices() three times faster by some hardwired code.
# 2006-11-28
# o Added readUnits() to override caching mechanism of superclasses.
# 2006-09-20
# o BUG FIX: Typo. Remove an argument but tried to use inside.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-12
# o Updated and probably working. When combining alleles, the names of the
#   groups returned consist of the allele A and allele group names.
# 2006-09-11
# o Created.
############################################################################
