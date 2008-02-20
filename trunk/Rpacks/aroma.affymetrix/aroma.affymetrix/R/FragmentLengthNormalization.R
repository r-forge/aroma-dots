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
#   \item{targetFunctions}{An optional list of @functions.  
#     For each enzyme there is one target function to which all arrays
#     should be normalized to.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
#   \item{shift}{An optional amount the data points should be shifted
#      (translated).}
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
setConstructorS3("FragmentLengthNormalization", function(dataSet=NULL, ..., targetFunctions=NULL, subsetToFit="-XY", shift=0) {
  extraTags <- NULL;

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

  # Argument 'targetFunctions':
  if (!is.null(targetFunctions)) {
    if (!is.list(targetFunctions)) {
      throw("Argument 'targetFunctions' is not a list: ", class(targetFunctions)[1]);
    }
    
    # Validate each element
    for (kk in seq(along=targetFunctions)) {
      if (!is.function(targetFunctions[[kk]])) {
        throw("One element in 'targetFunctions' is not a function: ", class(targetFunctions[[kk]])[1]);
      }
    }
  }

  # Argument 'subsetToFit':
  if (is.null(subsetToFit)) {
  } else if (is.character(subsetToFit)) {
    if (subsetToFit %in% c("-X", "-Y", "-XY")) {
    } else {
      throw("Unknown value of argument 'subsetToFit': ", subsetToFit);
    }
    extraTags <- c(extraTags, subsetToFit=subsetToFit);
  } else {
    cdf <- getCdf(dataSet);
    subsetToFit <- Arguments$getIndices(subsetToFit, 
                                        range=c(1, nbrOfUnits(cdf)));
    subsetToFit <- unique(subsetToFit);
    subsetToFit <- sort(subsetToFit);
  }

  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));


  extend(ChipEffectTransform(dataSet, ...), "FragmentLengthNormalization", 
    .subsetToFit = subsetToFit,
    .targetFunctions = targetFunctions,
    .extraTags = extraTags,
    shift = shift
  )
})



setMethodS3("getAsteriskTags", "FragmentLengthNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Add class-specific tags
  shift <- as.integer(round(this$shift));
  if (shift != 0) {
    tags <- c(tags, sprintf("%+d", shift));
  }

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, private=TRUE)


setMethodS3("clearCache", "FragmentLengthNormalization", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".targetFunctions")) {
    this[[ff]] <- NULL;
  }

  # Then for this object 
  NextMethod("clearCache", object=this, ...);
})


setMethodS3("getParameters", "FragmentLengthNormalization", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, expand=expand, ...);

  # Get parameters of this class
  params <- c(params, list(
    subsetToFit = this$.subsetToFit,
    .targetFunctions = this$.targetFunctions,
    shift <- this$shift
  ));


  # Expand?
  if (expand) {
    subsetToFit <- getSubsetToFit(this);
  }

  params;
}, private=TRUE)


setMethodS3("getCdf", "FragmentLengthNormalization", function(this, ...) {
  inputDataSet <- getInputDataSet(this);
  getCdf(inputDataSet);
})


setMethodS3("getOutputDataSet", "FragmentLengthNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting input data set");
  ces <- getInputDataSet(this);
  verbose && cat(verbose, "Class: ", class(ces)[1]);
  verbose && exit(verbose);

  verbose && enter(verbose, "Getting output data set for ", class(this)[1]);

  args <- list(generic="getOutputDataSet", object=this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Inherit certain arguments from the input data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(ces, "CnChipEffectSet"))
    args$combineAlleles <- ces$combineAlleles;
  if (inherits(ces, "SnpChipEffectSet"))
    args$mergeStrands <- ces$mergeStrands;

  verbose && cat(verbose, "Calling NextMethod() with arguments:");
  verbose && str(verbose, args);

  args$verbose <- less(verbose, 10);
  res <- do.call("NextMethod", args);

  # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
  if (inherits(res, "SnpChipEffectSet")) {
    verbose && enter(verbose, "Carrying down parameters for ", class(res)[1]);

    res$mergeStrands <- ces$mergeStrands;
    if (inherits(res, "CnChipEffectSet")) {
      res$combineAlleles <- ces$combineAlleles;
    }
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  res;
})

setMethodS3("getSubsetToFit", "FragmentLengthNormalization", function(this, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Cached?
  units <- this$.units;
  if (!is.null(units) && !force)
    return(units);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying all potential units, i.e. SNPs and CN probes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units that are SNP and CN probes");

  # Get SNP information
  cdf <- getCdf(this);
  si <- getSnpInformation(cdf);

  # Identify all SNP and CN units (==potential units to be fitted)
  verbose && enter(verbose, "Identifying SNPs and CN probes");
  units <- indexOf(cdf, "^(SNP|CN)");
  verbose && str(verbose, units);
  verbose && exit(verbose);

  # Keep only those for which we have PCR fragmenth-length information
  # for at least one enzyme
  verbose && enter(verbose, "Reading fragment lengths");
  fl <- getFragmentLengths(si, units=units);
  keep <- rep(FALSE, nrow(fl));
  for (ee in seq(length=ncol(fl))) {
    keep <- keep | is.finite(fl[,ee]);
  }
  units <- units[keep];
  verbose && printf(verbose, "Number of SNP/CN units without fragment-length details: %d out of %d (%.1f%%)\n", sum(!keep), length(keep), 100*sum(!keep)/length(keep));

  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset with a prespecified set of units?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subsetToFit <- this$.subsetToFit;
  if (is.character(subsetToFit)) {
    if (subsetToFit %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Identify subset of units from genome information");
      verbose && cat(verbose, "subsetToFit: ", subsetToFit);

      # Look up in cache
      subset <- this$.subsetToFitExpanded;
      if (is.null(subset)) {
        # Get the genome information (throws an exception if missing)
        gi <- getGenomeInformation(cdf);
        verbose && print(verbose, gi);
  
        # Identify units to be excluded
        if (subsetToFit == "-X") {
          subset <- getUnitsOnChromosome(gi, 23, .checkArgs=FALSE);
        } else if (subsetToFit == "-Y") {
          subset <- getUnitsOnChromosome(gi, 24, .checkArgs=FALSE);
        } else if (subsetToFit == "-XY") {
          subset <- getUnitsOnChromosome(gi, 23:24, .checkArgs=FALSE);
        }
  
        verbose && cat(verbose, "Units to exclude: ");
        verbose && str(verbose, subset);
  
        # The units to keep
        subset <- setdiff(1:nbrOfUnits(cdf), subset);
  
        verbose && cat(verbose, "Units to include: ");
        verbose && str(verbose, subset);

        # Store
        this$.subsetToFitExpanded <- subset;
      }

      subsetToFit <- subset;
      rm(subset);

      verbose && exit(verbose);
    }
  }

  if (!is.null(subsetToFit)) {
    # A fraction subset?
    if (length(subsetToFit) == 1 && 0 < subsetToFit && subsetToFit < 1) {
      keep <- seq(from=1, to=length(units), length=subsetToFit*length(units));
    } else {
      keep <- which(units %in% subsetToFit);
    }

    # Make sure to keep data points at the tails too
    extremeUnits <- c();
    for (ee in seq(length=ncol(fl))) {
      extremeUnits <- c(extremeUnits, which.min(fl[,ee]), which.max(fl[,ee]));
    }
    keep <- c(keep, extremeUnits);
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

  verbose && exit(verbose);

  units;
}, private=TRUE)



setMethodS3("getTargetFunction", "FragmentLengthNormalization", function(this, ...) {
  getTargetFunctions(this, ...);
}, deprecated=TRUE)

setMethodS3("getTargetFunctions", "FragmentLengthNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  fcns <- this$.targetFunctions;
  if (is.null(fcns) || force) {
    verbose && enter(verbose, "Estimating target prediction function");

    # Get the SNP information annotation
    cdf <- getCdf(this);
    si <- getSnpInformation(cdf);
    rm(cdf); # Not needed anymore

    # Get target set
    ces <- getInputDataSet(this);
    verbose && enter(verbose, "Get average signal across arrays");
    ceR <- getAverageFile(ces, force=force, verbose=less(verbose));
    rm(ces); # Not needed anymore
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Get units to fit
    units <- getSubsetToFit(this);

    # Get target log2 signals for SNPs
    data <- getDataFlat(ceR, units=units, fields="theta", verbose=less(verbose));
    rm(ceR); # Not needed anymore
    units <- data[,"unit"];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    yR <- data[,"theta"];
    rm(data); # Not needed anymore

    # Shift?
    shift <- this$shift;
    if (shift != 0) {
      yR <- yR + shift;
    }
    
    yR <- log2(yR);
    verbose && cat(verbose, "Signals:");
    verbose && str(verbose, yR);
    
    # Get PCR fragment lengths for these
    fl <- getFragmentLengths(si, units=units);
    rm(si, units); # Not needed anymore
    verbose && cat(verbose, "Fragment lengths:");
    verbose && str(verbose, fl);

    # Fit lowess function
    verbose && enter(verbose, "Fitting target prediction function to each enzyme exclusively");
    okYR <- is.finite(yR);
    hasFL <- is.finite(fl);
    nbrOfEnzymes <- ncol(fl);
    allEnzymes <- seq(length=nbrOfEnzymes);

    fits <- list();
    for (ee in allEnzymes) {
      verbose && enter(verbose, "Enzyme #", ee, " of ", nbrOfEnzymes);

      # Fit only to units with known length and non-missing data points.
      ok <- (hasFL[,ee] & okYR);

      # Exclude multi-enzyme units
      for (ff in setdiff(allEnzymes, ee))
        ok <- ok & !hasFL[,ff];

      # Sanity check
      if (sum(ok) == 0) {
        throw("Cannot fit target function to enzyme, because there are no (finite) data points that are unique to this enzyme: ", ee);
      }

      # Fit fragment-length effect to single-enzyme units
      fit <- lowess(fl[ok,ee], yR[ok]);
      class(fit) <- "lowess";

      rm(ok);

      fits[[ee]] <- fit;

      rm(fit);

      verbose && exit(verbose);
    }

    # Remove as many promises as possible
    rm(fcns, nbrOfEnzymes, allEnzymes, fl, yR, okYR, hasFL);

    # Create a target prediction function for each enzyme
    fcns <- vector("list", length(fits));
    for (ee in seq(along=fits)) {
      fcns[[ee]] <- function(x, ...) {
        predict(fits[[ee]], x, ...);  # Dispatched predict.lowess().
      }
    }
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);

    this$.targetFunctions <- fcns;
  }

  fcns;
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
    verbose && enter(verbose, "Getting output data set");
    outputSet <- getOutputDataSet(this, verbose=less(verbose, 10));
    verbose && exit(verbose);
    verbose && exit(verbose);
    return(invisible(outputSet));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ces <- getInputDataSet(this);

  # Get SNP units
  cdf <- getCdf(ces);
  subsetToUpdate <- indexOf(cdf, "^(SNP|CN)");

  verbose && enter(verbose, "Retrieving SNP information annotations");
  si <- getSnpInformation(cdf);
  verbose && cat(verbose, "Number of enzymes: ", nbrOfEnzymes(si));
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying the subset used to fit normalization function(s)");
  # Get subset to fit
  subsetToFit <- getSubsetToFit(this, verbose=less(verbose));
  verbose && str(verbose, subsetToFit);
  verbose && exit(verbose);

  shift <- this$shift;
  verbose && cat(verbose, "Shift: ", shift);

  # Get (and create) the output path
  path <- getPath(this);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fl <- NULL;
  targetFcns <- NULL;
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

    # Get unit-to-cell? (for optimized reading)
    if (is.null(map)) {
      # Only loaded if really needed.
      verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
      map <- getCellMap(ce, units=subsetToUpdate, verbose=less(verbose));
      verbose && str(verbose, map);
      verbose && exit(verbose);
    }

    if (is.null(fl)) {
      # For the subset to be fitted, get PCR fragment lengths (for all enzymes) 
      fl <- getFragmentLengths(si, units=map[,"unit"]);

      # Get the index in the data vector of subset to be fitted.
      # Note: match() only returns first match, which is why we do
      # it this way.
      subset <- match(map[,"unit"], subsetToFit);
      subset <- subset[!is.na(subset)];
      subset <- match(subsetToFit[subset], map[,"unit"]);
    }

    if (is.null(targetFcns)) {
      # Only loaded if really needed.
      # Retrieve/calculate the target function
      targetFcns <- getTargetFunctions(this, verbose=less(verbose));
    }

    # Get target log2 signals for all SNPs to be updated
    verbose && enter(verbose, "Getting signals");
    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
    verbose && exit(verbose);


    # Extract the values to fit the normalization function
    verbose && enter(verbose, "Normalizing log2 signals");
    y <- data[,"theta"];

    # Shift?
    if (shift != 0)
      y <- y + shift;

    # Fit on the log2 scale
    y <- log2(y);

    verbose && cat(verbose, "Log2 signals:");
    verbose && str(verbose, y);
    y <- normalizeFragmentLength(y, fragmentLengths=fl, 
                             targetFcns=targetFcns, subsetToFit=subset, ...);
    
    # Store results on the intensity scale
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
# 2008-02-18
# o Added 'shift' to FragmentLengthNormalization, cf ProbeLevelModel.
# 2007-12-01
# o Added getAsteriskTag() to FragmentLengthNormalization.
# o Similar to AllelicCrosstalkCalibration, the constructor argument
#   'subsetToFit' of FragmentLengthNormalization accept "-XY" (and "-X" and
#   "-Y") to specify the set of units to fit the model over to be all units
#   that are not on ChrX or ChrY.
# 2007-11-19
# o Updated getSubsetToFit() to handle chip types with multiple enzymes.
# o Updated methods to give an error if chip types with more than one
#   enzyme is tried to be normalized.
# 2007-09-16
# o Added clearCache() to FragmentLengthNormalization such that cached
#   target distributions can be cleared.
# 2007-09-12
# o BUG FIX: getSubsetToFit() of FragmentLengthNormalization would only
#   return SNP units, but not CN units which are available on the newer
#   chip types.  Similarly, process() would only update SNPs, but not
#   CN units.
# o Now getOutputDataSet() of FragmentLengthNormalization set and pass down
#   'mergeStrands' and 'combineAlleles' to ditto of the super class, if
#   applicable.  This way we avoid having to infer those arguments from
#   the contents of the files.
# 2007-02-20
# o Now FragmentLengthNormalization should handle cases with more than one
#   chip effect per unit, e.g. when mergeStrands=FALSE.
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
