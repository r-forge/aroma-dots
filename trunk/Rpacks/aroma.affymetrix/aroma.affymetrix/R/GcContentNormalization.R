###########################################################################/**
# @RdocClass GcContentNormalization
#
# @title "The GcContentNormalization class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{A @see "CnChipEffectSet".}
#   \item{...}{Additional arguments passed to the constructor of 
#     @see "ChipEffectTransform".}
#   \item{targetFunction}{A @function. The target function to which all arrays
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
#   This class requires an Affymetrix TSV file for the chip type to be
#   normalized.
# }
#
# \examples{\dontrun{
# }}
#
# @author
#*/###########################################################################
setConstructorS3("GcContentNormalization", function(dataSet=NULL, ..., targetFunction=NULL, subsetToFit=NULL) {
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

#    if (dataSet$mergeStrands != TRUE) {
#      throw("Currently only non-strands specific copy-number chip effects can be normalized, i.e. 'mergeStrands' must be TRUE");
#    }
  }

  if (!is.null(targetFunction)) {
    if (!is.function(targetFunction)) {
      throw("Argument 'targetFunction' is not a function: ", class(targetFunction)[1]);
    }
  }

  extend(ChipEffectTransform(dataSet, ...), "GcContentNormalization", 
    .subsetToFit = subsetToFit,
    .targetFunction = targetFunction
  )
})



setMethodS3("getParameters", "GcContentNormalization", function(this, ...) {
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


setMethodS3("getCdf", "GcContentNormalization", function(this, ...) {
  inputDataSet <- getInputDataSet(this);
  getCdf(inputDataSet);
})


setMethodS3("getOutputDataSet", "GcContentNormalization", function(this, ...) {
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

setMethodS3("getSubsetToFit", "GcContentNormalization", function(this, force=FALSE, ...) {
  # Cached?
  units <- this$.units;
  if (!is.null(units) && !force)
    return(units);

  # Identify all SNP (& CN) units
  cdf <- getCdf(this);
#  units <- indexOf(cdf, "SNP_");
  types <- getUnitTypes(cdf, ...);
  units <- which(types == 2 | types == 5 | types == 8);

  # Keep only those for which we have GC contents information
  tsv <- AffymetrixTsvFile$fromChipType(getChipType(cdf, fullname=FALSE));
  gcContents <- getGc(tsv, units=units);
  keep <- is.finite(gcContents);
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
    keep <- c(keep, which.min(gcContents), which.max(gcContents));
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



setMethodS3("getTargetFunction", "GcContentNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
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

    # Get the Affymetrix TSV file
    cdf <- getCdf(this);
    tsv <- AffymetrixTsvFile$fromChipType(getChipType(cdf, fullname=FALSE));

    # Get target set
    ces <- getInputDataSet(this);
    verbose && enter(verbose, "Get average signal across arrays");
    ceR <- getAverageFile(ces, force=force, verbose=less(verbose));
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Get units to fit
    units <- getSubsetToFit(this);

    # Get target log2 signals for SNPs
    data <- getDataFlat(ceR, units=units, fields="theta", verbose=less(verbose));
    units <- data[,"unit"];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    yR <- data[,"theta"];
    rm(data); # Not needed anymore
    yR <- log2(yR);
    verbose && cat(verbose, "Signals:");
    verbose && str(verbose, yR);
    
    # Get PCR fragment lengths for these
    gcContents <- getGc(tsv, units=units);
    verbose && cat(verbose, "GC content:");
    verbose && str(verbose, gcContents);

    # Fit lowess function
    verbose && enter(verbose, "Fitting target prediction function");
    ok <- (is.finite(gcContents) & is.finite(yR));
    fit <- lowess(gcContents[ok], yR[ok]);
    class(fit) <- "lowess";

    # Remove as many promises as possible
    rm(fcn, ces, ceR, units, cdf, tsv, gc, yR, ok);

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
setMethodS3("process", "GcContentNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
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

  # Get SNP (& CN) units
  cdf <- getCdf(ces);
#  subsetToUpdate <- indexOf(cdf, "SNP_");
  types <- getUnitTypes(cdf, ...);
  subsetToUpdate <- which(types == 2 | types == 5 | types == 8);

  verbose && enter(verbose, "Retrieving Affymetrix TSV file");
  tsv <- AffymetrixTsvFile$fromChipType(getChipType(cdf, fullname=FALSE));
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying the subset used to fit normalization function");
  # Get subset to fit
  subsetToFit <- getSubsetToFit(this, verbose=less(verbose));
  verbose && str(verbose, subsetToFit);
  verbose && exit(verbose);

  # Get (and create) the output path
  path <- getPath(this);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  gcContents <- NULL;
  targetFcn <- NULL;
  map <- NULL;
  nbrOfArrays <- nbrOfArrays(ces);
  res <- vector("list", nbrOfArrays);
  for (kk in seq_len(nbrOfArrays)) {
    ce <- getFile(ces, kk);
    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')",
                                            kk, nbrOfArrays, getName(ce)));

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

    if (is.null(gcContents)) {
      # Get PCR fragment lengths for the subset to be fitted
      gcContents <- getGc(tsv, units=map[,"unit"]);

      # Get the index in the data vector of subset to be fitted.
      # Note: match() only returns first match, which is why we do
      # it this way.
      subset <- match(map[,"unit"], subsetToFit);
      subset <- subset[!is.na(subset)];
      subset <- match(subsetToFit[subset], map[,"unit"]);
    }

    if (is.null(targetFcn)) {
      # Only loaded if really needed.
      # Retrieve/calculate the target function
      targetFcn <- getTargetFunction(this, verbose=less(verbose));
    }

    # Get target log2 signals for all SNPs to be updated
    verbose && enter(verbose, "Getting signals");
    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
    verbose && exit(verbose);


    # Extract the values to fit the normalization function
    verbose && enter(verbose, "Normalizing log2 signals");
    y <- log2(data[,"theta"]);
    y <- normalizeFragmentLength(y, fragmentLengths=gcContents, 
                             targetFcn=targetFcn, subsetToFit=subset, ...);
    y <- 2^y;
    verbose && exit(verbose);

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");
    ceN <- createFrom(ce, filename=pathname, path=NULL, verbose=less(verbose));
    verbose && exit(verbose);

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
    data[,"theta"] <- y;
    rm(y);
    updateDataFlat(ceN, data=data, verbose=less(verbose));
    rm(data);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    res[[kk]] <- ceN;

    verbose && exit(verbose);
  } # for (kk in ...)

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
# 2008-02-21
# o Now getSubsetToFit() and process() not only processes SNPs but also
#   CN probes.
# 2007-04-02
# o Created from FragmentLengthNormalization.R.
############################################################################
