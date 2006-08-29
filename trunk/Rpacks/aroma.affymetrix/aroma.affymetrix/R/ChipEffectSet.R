###########################################################################/**
# @RdocClass ChipEffectSet
#
# @title "The ChipEffectSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ParameterCelFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getChipEffects()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("ChipEffectSet", function(..., model=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  extend(AffymetrixCelSet(...), "ChipEffectSet",
    "cached:.firstCells" = NULL,
    model = model
  )
})


setMethodS3("getFirstCellIndices", "ChipEffectSet", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose': 
  verbose <- Arguments$getVerbose(verbose);

  res <- this$.firstCells;
  if (is.null(res)) {
    stratifyBy <- switch(this$model, pm="pm");
    cdf <- getCdf(this);
    res <- getFirstCellIndices(cdf, units=NULL, ..., stratifyBy=stratifyBy, verbose=verbose);
    this$.firstCells <- res;
  }

  # Subset?
  if (!is.null(units))
    res <- res[units];

  res;
}, protected=TRUE)



setMethodS3("readUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, ...) {
  if (is.null(cdf)) {
    # Use only the first cell in each unit group.
    cdf <- getFirstCellIndices(this, units=units, ...);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  stratifyBy <- switch(this$model, pm="pm");
  res <- NextMethod("readUnits", this, units=cdf, ..., stratifyBy=stratifyBy);

  # Get first chip-effect file and use that to decode the read structure
  ce <- as.list(this)[[1]];
  res <- decode(ce, res);

  res;
});


setMethodS3("updateUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # Argument 'verbose': 
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF structure for all chip-effect files
  if (is.null(cdf)) {
    # Use only the first cell in each unit group.
    cdf <- getFirstCellIndices(this, units=units, ...);
  }

  # Update each file one by one
  n <- length(this);
  names <- getNames(this);
  for (kk in seq(this)) {
    verbose && enter(verbose, sprintf("Array #%d of %d: %s", kk, n, names[kk]));
    ce <- as.list(this)[[kk]];

    verbose && enter(verbose, "Extracting estimates");
    dataOne <- lapply(data, function(groups) {
      lapply(groups, function(group) {
        # theta = group$theta[kk] = ...
        # stdvs = 1 (default for now)
        list(theta=.subset(.subset2(group, "theta"), kk), stdvs=1);
      })
    })
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating file");
    updateUnits(ce, cdf=cdf, data=dataOne);
    verbose && exit(verbose);
    verbose && exit(verbose);
  }
}, protected=TRUE);


setMethodS3("getAverageFile", "ChipEffectSet", function(this, indices="remaining", cellsPerChunk=100, ...) {
  # Argument 'indices':
  if (identical(indices, "remaining")) {
  } else if (is.null(indices)) {
    # Update only cells which stores values
    indices <- getFirstCellIndices(this);
    indices <- unlist(indices, use.names=FALSE);
  }

  NextMethod("getAverageFile", this, indices=indices, cellsPerChunk=cellsPerChunk, ...);
})


############################################################################
# HISTORY:
# 2006-08-28
# o Added getAverageFile() so that only cells that store actual chip-effect
#   estimates are averaged.
# 2006-08-26
# o Created.
############################################################################
