###########################################################################/**
# @RdocClass FragmentLengthNormalization
#
# @title "The FragmentLengthNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for PCR 
#  fragment length effects on total copy-number chip-effect estimates.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{A @see "CnChipEffectSet".}
#   \item{...}{Additional arguments passed to the constructor of 
#     @see "ChipEffectTransform".}
#   \item{targetFunction}{A @function.  The target function to which all arrays
#     should be normalized to.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \section{Requirements}{
#   This class requires a SNP information annotation file for the 
#   chip type to be normalized.
# }
#
# \examples{\dontrun{
# }}
#
# @author
#*/###########################################################################
setConstructorS3("FragmentLengthNormalization", function(dataSet=NULL, ..., targetFunction=NULL, subsetToFit=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "CnChipEffectSet"))
      throw("Argument 'dataSet' is not an CnChipEffectSet object: ", class(dataSet));

    if (dataSet$combineAlleles != TRUE) {
      throw("Currently only total copy-number chip effects can be normalized, i.e. 'combineAlleles' must be TRUE");
    }

    if (dataSet$mergeStrands != TRUE) {
      throw("Currently only non-strands specific copy-number chip effects can be normalized, i.e. 'mergeStrands' must be TRUE");
    }
  }

  if (!is.null(targetFunction)) {
    if (!is.function(targetFunction)) {
      throw("Argument 'targetFunction' is not a function: ", class(targetFunction)[1]);
    }
  }

  extend(ChipEffectTransform(dataSet, ...), "FragmentLengthNormalization", 
    .subsetToFit = subsetToFit,
    .targetFunction = targetFunction
  )
})



setMethodS3("getParameters", "FragmentLengthNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  # Get parameters of this class
  params2 <- list(
    subsetToFit = this$.subsetToFit,
    .targetFunction = this$.targetFunction
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)


setMethodS3("getCdf", "FragmentLengthNormalization", function(this, ...) {
  inputDataSet <- getInputDataSet(this);
  getCdf(inputDataSet);
})


setMethodS3("getOutputDataSet", "FragmentLengthNormalization", function(this, ...) {
  res <- NextMethod(generic="getOutputDataSet", object=this, ...);

  # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
  if (inherits(res, "SnpChipEffectSet")) {
    ces <- getInputDataSet(this);
    res$mergeStrands <- ces$mergeStrands;
    if (inherits(res, "CnChipEffectSet")) {
      res$combineAlleles <- ces$combineAlleles;
    }
  }

  res;
})

setMethodS3("getSubsetToFit", "FragmentLengthNormalization", function(this, force=FALSE, ...) {
  # Cached?
  units <- this$.units;
  if (!is.null(units) && !force)
    return(units);

  # Identify all SNP units
  cdf <- getCdf(this);
  units <- indexOf(cdf, "^SNP");

  # Keep only those for which we have PCR fragmenth-length information
  si <- getSnpInformation(cdf);
  fl <- getFragmentLengths(si, units=units);
  keep <- is.finite(fl);
  units <- units[keep];

  # Fit to a subset of the units?
  subsetToFit <- this$.subsetToFit;
  if (!is.null(subsetToFit)) {
    # A fraction subset?
    if (length(subsetToFit) == 1 && 0 < subsetToFit && subsetToFit < 1) {
      keep <- seq(from=1, to=length(units), length=subsetToFit*length(units));
    } else {
      keep <- which(units %in% subsetToFit);
    }

    # Make sure to keep data points at the tails too
    keep <- c(keep, which.min(fl), which.max(fl));
    keep <- unique(keep);

    # Now filter
    units <- units[keep];
    rm(keep);
  }

  # Sort units
  units <- sort(units);

  # Assert correctness
  units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));

  # Cache
  this$.units <- units;

  units;
}, private=TRUE)



setMethodS3("getTargetFunction", "FragmentLengthNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  fcn <- this$.targetFunction;
  if (is.null(fcn) || force) {
    verbose && enter(verbose, "Estimating target prediction function");

    # Get target set
    ces <- getInputDataSet(this);
    verbose && enter(verbose, "Get average signal across arrays");
    ceR <- getAverageFile(ces, indices=NULL, force=force, verbose=less(verbose));
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Get units to fit
    units <- getSubsetToFit(this);

    # Get PCR fragment lengths for these
    cdf <- getCdf(this);
    si <- getSnpInformation(cdf);
    fl <- getFragmentLengths(si, units=units);

    # Get target log2 signals for SNPs
    yR <- getDataFlat(ceR, units=units, fields="theta", verbose=less(verbose));
    yR <- yR[,"theta"];
    yR <- log2(yR);

    # Fit lowess function
    verbose && enter(verbose, "Fitting target prediction function");
    ok <- (is.finite(fl) & is.finite(yR));
    fit <- lowess(fl[ok], yR[ok]);
    class(fit) <- "lowess";

    # Remove as many promises as possible
    rm(fcn, ces, ceR, units, cdf, si, fl, yR, ok);

    # Create target prediction function
    fcn <- function(x, ...) {
      predict(fit, x, ...);  # Dispatched predict.lowess().
    }
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);

    this$.targetFunction <- fcn;
  }

  fcn;
}, private=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized, 
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @examples "../incl/normalizeQuantile.Rex"
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "FragmentLengthNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing set for PCR fragment-length effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputSet <- getOutputDataSet(this);
    return(invisible(outputSet));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ces <- getInputDataSet(this);

  # Get SNP units
  cdf <- getCdf(ces);
  subsetToUpdate <- indexOf(cdf, "^SNP");

  verbose && enter(verbose, "Retrieving SNP information annotations");
  cdf <- getCdf(this);
  si <- getSnpInformation(cdf);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying the subset to be fitted");
  # Get subset to fit
  subsetToFit <- getSubsetToFit(this, verbose=less(verbose));
  verbose && exit(verbose);

  # Get (and create) the output path
  path <- getPath(this);
  mkdirs(path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fl <- NULL;
  targetFcn <- NULL;
  map <- NULL;
  res <- vector("list", nbrOfArrays(ces));
  for (kk in seq(length=nbrOfArrays(ces))) {
    ce <- getFile(ces, kk);
    verbose && enter(verbose, sprintf("Array #%d (%s)", kk, getName(ce)));

    filename <- getFilename(ce);
    pathname <- filePath(path, filename);
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized. Skipping.");
      ceN <- fromFile(ce, pathname);

      # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
      if (inherits(ce, "SnpChipEffectFile")) {
        ceN$mergeStrands <- ce$mergeStrands;
        if (inherits(ce, "CnChipEffectFile")) {
          ceN$combineAlleles <- ce$combineAlleles;
        }
      }

      # CDF inheritance
      setCdf(ceN, cdf);

      res[[kk]] <- ceN;
      verbose && exit(verbose);
      next;
    }

    # Get unit-to-cell (for optimized reading)?
    if (is.null(map)) {
      # Only loaded if really needed.
      verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
      map <- getCellMap(ce, units=subsetToUpdate, verbose=less(verbose));
      verbose && str(verbose, map);
      verbose && exit(verbose);
    }

    if (is.null(targetFcn)) {
      # Only loaded if really needed.
      # Retrieve/calculate the target function
      targetFcn <- getTargetFunction(this, verbose=less(verbose));
    }

    if (is.null(fl)) {
      # Get PCR fragment lengths for the subset to be fitted
      fl <- getFragmentLengths(si, units=subsetToUpdate);
    }

    # Get target log2 signals for all SNPs to be updated
    verbose && enter(verbose, "Getting signals");
    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
    verbose && exit(verbose);

    # Extract the values to fit the normalization function
    verbose && enter(verbose, "Normalizing log2 signals");
    y <- data[,"theta"];
    subset <- match(subsetToFit, subsetToUpdate);
    yN <- normalizeFragmentLength(log2(y), fragmentLengths=fl, 
                             targetFcn=targetFcn, subsetToFit=subset, ...);
    rm(y);
    yN <- 2^yN;
    verbose && exit(verbose);

    # Copy CEL file and update the copy
    verbose && enter(verbose, "Copying source CEL file");
    copyCel(from=getPathname(ce), to=pathname);
    verbose && exit(verbose);

    # Defining normalized object
    ceN <- fromFile(ce, pathname);

    # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
    if (inherits(ce, "SnpChipEffectFile")) {
      ceN$mergeStrands <- ce$mergeStrands;
      if (inherits(ce, "CnChipEffectFile")) {
        ceN$combineAlleles <- ce$combineAlleles;
      }
    }

    # CDF inheritance
    setCdf(ceN, cdf);

    verbose && enter(verbose, "Storing normalized signals");
    data[,"theta"] <- yN;
    rm(yN);
    updateDataFlat(ceN, data=data, verbose=less(verbose));
    rm(data);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    res[[kk]] <- ceN;

    verbose && exit(verbose);
  }

  # Create the output set (ad hoc for now so that we keep parameter too)
  outputSet <- clone(ces);
  outputSet$files <- res;
  clearCache(outputSet);

  # Update the output data set
  this$outputSet <- outputSet;

  verbose && exit(verbose);
  
  outputSet;
})

############################################################################
# HISTORY:
# 2007-01-16
# o BUG FIX: Forgot to clear the cache after cloning data set in process().
#   This would cause getAverage() to return a cached averaged from the
#   non-normalized data set.
# 2007-01-07
# o Now chip-effect parameters are carried over to the output set too.
# o BUG FIX: process(): Forgot to skip to next array in for loop if an 
#   array was detected to be already normalized. Generated a "file already
#   exists" error.
# o Added garbage collection after each array been normalized.
# 2007-01-04
# o BUG FIX: process() gave an error if the data set was already done.
# 2006-12-08
# o Now this class inherits from the ChipEffectPreprocessing class.
# o Now this pre-processor output results to plmData/.
# 2006-11-28
# o Created from QuantileNormalizer.R.
############################################################################
