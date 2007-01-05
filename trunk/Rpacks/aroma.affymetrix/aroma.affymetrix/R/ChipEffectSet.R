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
#   \item{model}{The specific type of model, e.g. \code{"pm"}.}
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


setMethodS3("getChipEffectFileClass", "ChipEffectSet", function(static, ...) {
  ChipEffectFile;
}, static=TRUE, private=TRUE)


setMethodS3("fromFiles", "ChipEffectSet", function(static, ..., pattern=",chipEffects[.](c|C)(e|E)(l|L)$", fileClass=NULL) {
  # Argument 'fileClass':
  if (is.null(fileClass))
    fileClass <- gsub("Set$", "File", class(static)[1]);

  fromFiles.AffymetrixFileSet(static, ..., pattern=pattern, fileClass=fileClass);
}, static=TRUE)


setMethodS3("fromDataSet", "ChipEffectSet", function(static, dataSet, path, name=getName(dataSet), ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the ChipEffectFile class specific for this set
  clazz <- getChipEffectFileClass(static);

  verbose && enter(verbose, "Retrieving chip-effects from data set");
  ces <- vector("list", length(dataSet));
  verbose && cat(verbose, "Data set: ", name);
  cdf <- NULL;
  for (kk in seq(dataSet)) {
    df <- getFile(dataSet, kk);
    verbose && enter(verbose, 
                           sprintf("Retrieving chip-effect #%d of %d (%s)",
                                               kk, length(ces), getName(df)));
    ce <- clazz$fromDataFile(df, path=path, name=name, cdf=cdf, ..., 
                                                       verbose=less(verbose));
    if (is.null(cdf)) {
      verbose && enter(verbose, "Requiring the CDF for the chip-effect file");
      cdf <- getCdf(ce);
      verbose && exit(verbose);
    }
    ces[[kk]] <- ce;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Create an ChipEffectSet
  newInstance(static, ces);
})


setMethodS3("getCellIndices", "ChipEffectSet", function(this, ...) {
  # Use the first chip-effect file to get the CDF structure.
  # Note: Ideally we want to define a special CDF class doing this
  # instead of letting the data file do this. /HB 2006-12-18
  ce <- getFile(this, 1);
  getCellIndices(ce, ...);
})


setMethodS3("readUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading chip effects unit by unit for ", nbrOfArrays(this), " arrays");

  if (is.null(cdf)) {
    verbose && enter(verbose, "Getting cell indices from CDF");
    cdf <- getCellIndices(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  verbose && enter(verbose, "Calling readUnits() in superclass");
  res <- NextMethod("readUnits", this, units=cdf, ..., verbose=less(verbose));
  verbose && exit(verbose);

  # Get first chip-effect file and use that to decode the read structure
  # This takes some time for a large number of units /HB 2006-10-04
  ce <- getFile(this, 1);
  res <- decode(ce, res, verbose=less(verbose));

  verbose && exit(verbose);

  res;
})


setMethodS3("updateUnits", "ChipEffectSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # Argument 'verbose': 
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF structure for all chip-effect files
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Update each file one by one
  n <- length(this);
  verbose && enter(verbose, "Updating ", n, " chip-effect files");
  names <- getNames(this);
  verbose <- less(verbose);
  for (kk in seq(this)) {
    verbose && enter(verbose, sprintf("Array #%d of %d: %s", kk, n, names[kk]));
    ce <- as.list(this)[[kk]];

    verbose <- less(verbose, 50);
    verbose && enter(verbose, "Extracting estimates");  # 3-4s
    dataOne <- lapply(data, FUN=lapply, function(group) {
      # theta = group$theta[kk] = ...
      # stdvs = group$sdTheta[kk] = ...
      list(
        theta=.subset(.subset2(group, "theta"), kk), 
        sdTheta=.subset(.subset2(group, "sdTheta"), kk),
        thetaOutliers=.subset(.subset2(group, "thetaOutliers"), kk)
      );
    });
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating file");  # 6-7s ~98% in encode()
    updateUnits(ce, cdf=cdf, data=dataOne, verbose=less(verbose, 50));
    verbose && exit(verbose);
    verbose <- more(verbose, 50);

    verbose && exit(verbose);
  } # for (kk ...)
  verbose <- more(verbose);
  verbose && exit(verbose);
}, protected=TRUE)



setMethodS3("getAverageFile", "ChipEffectSet", function(this, ..., indices="remaining") {
  # Argument 'indices':
  if (identical(indices, "remaining")) {
  } else if (is.null(indices)) {
    # Update only cells which stores values
    indices <- getCellIndices(this);
    indices <- unlist(indices, use.names=FALSE);
  }

  NextMethod(generic="getAverageFile", object=this, ..., indices=indices);
})


setMethodS3("findUnitsTodo", "ChipEffectSet", function(this, ...) {
  # Look into the last chip-effect file since that is updated last
  ce <- getFile(this, length(this));
  findUnitsTodo(ce, ...);
})

############################################################################
# HISTORY:
# 2006-11-22
# o Updated fromFiles() so it automagically finds the file class.
# 2006-10-??
# o Updated fromFiles() to search for filename with tags 'chipEffects'.
#   Before files with suffix "-chipEffects' was searched for.
# 2006-09-10
# o Added findUnitsTodo().
# o Starting to make use of specially design CDFs and CEL files for storing
#   chip effects.  This make getFirstCellIndices() obsolete.
# 2006-08-28
# o Added getAverageFile() so that only cells that store actual chip-effect
#   estimates are averaged.
# 2006-08-26
# o Created.
############################################################################
