setMethodS3("boxplotStats", "ChipEffectSet", function(this, type=c("NUSE", "RLE"), ...) {
  if (toupper(type) == "NUSE") {
    calculateNuseBoxplotStats(this, ...);    
  } else if (toupper(type) == "RLE") {
    calculateRleBoxplotStats(this, ...);    
  } else {
    calculateFieldBoxplotStats(this, field=type, ...);
  }
})




setMethodS3("calculateFieldBoxplotStats", "ChipEffectSet", function(this, field=c("theta", "sdTheta"), transform=NULL, arrays=NULL, subset=NULL, ..., merge=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfArrays <- nbrOfArrays(this);
  cdfMono <- getCdf(this);
  nbrOfUnits <- nbrOfUnits(cdfMono);

  # Argument 'field':
  field <- match.arg(field);

  # Argument 'transform':
  if (is.null(transform)) {
  } else if (is.function(transform)) {
  } else {
    throw("Argument 'transform' is not a function: ", class(transform)[1]);
  }


  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq(length=nbrOfArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, range=c(1, nbrOfArrays));
    nbrOfArrays <- length(arrays);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'subset':  
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, range=c(1, nbrOfUnits));
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }


  verbose && enter(verbose, "Calculating '", field, 
                  "' statistics for ", nbrOfArrays, " (specified) arrays");
  
  # For each file, calculate boxplot statistics
  stats <- list();
  for (kk in seq(along=arrays)) {
    array <- arrays[kk];
    cef <- getFile(this, array);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                          kk, getName(cef), nbrOfArrays));
    data <- extractMatrix(cef, units=units, fields=field);
    data <- as.vector(data);
    if (is.function(transform)) {
      data <- transform(data);
    }
    stats[[kk]] <- boxplot.stats(data);
    rm(data);
    verbose && exit(verbose);
  }
  names(stats) <- getNames(this)[arrays];
  verbose && exit(verbose);

  # Merge boxplot stats?
  if (merge) 
    stats <- mergeBoxplotStats(stats);

  attr(stats, "type") <- field;
  attr(stats, "transform") <- transform;

  stats;
}) # calculateFieldBoxplotStats()






# ... : additional arguments to bxp().
setMethodS3("calculateRleBoxplotStats", "ChipEffectSet", function(this, arrays=NULL, subset=NULL, ..., merge=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfArrays <- nbrOfArrays(this);
  cdfMono <- getCdf(this);
  nbrOfUnits <- nbrOfUnits(cdfMono);

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq(length=nbrOfArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, range=c(1, nbrOfArrays));
    nbrOfArrays <- length(arrays);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'subset':  
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, range=c(1, nbrOfUnits));
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }


  

  # get the vector of median stdvs
  verbose && enter(verbose, "Calculating average log chip effects");
  avg <- getAverageLog(this, field="intensities", indices=units, 
                                                         verbose=verbose);
  verbose && exit(verbose);

  medianLE <- getData(avg, indices=units, "intensities")$intensities;
  medianLE <- log2(medianLE);
  rm(avg);

  verbose && enter(verbose, "Calculating RLE statistics for ", nbrOfArrays, 
                                                    " (specified) arrays");
  
  # For each file, calculate boxplot statistics
  stats <- list();
  for (kk in seq(along=arrays)) {
    array <- arrays[kk];
    cef <- getFile(this, array);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                          kk, getName(cef), nbrOfArrays));

    theta <- getData(cef, indices=units, "intensities")$intensities;
    theta <- log2(theta);
    stats[[kk]] <- boxplot.stats(theta-medianLE);
    rm(theta);
    verbose && exit(verbose);
  }
  rm(medianLE, units);
  names(stats) <- getNames(this)[arrays];
  verbose && exit(verbose);

  # Merge boxplot stats?
  if (merge) 
    stats <- mergeBoxplotStats(stats);

  attr(stats, "type") <- "RLE";

  stats;
}) # calculateRleBoxplotStats()



setMethodS3("calculateNuseBoxplotStats", "ChipEffectSet", function(this, arrays=NULL, subset=NULL, ..., merge=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdfMono <- getCdf(this);
  nbrOfUnits <- nbrOfUnits(cdfMono);
  nbrOfArrays <- nbrOfArrays(this);

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq(length=nbrOfArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, range=c(1, nbrOfArrays));
    nbrOfArrays <- length(arrays);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'subset':  
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, range=c(1, nbrOfUnits));
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }

 

  # get the vector of median stdvs
  verbose && enter(verbose, "Extracting average standard errors across all arrays in the set");
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  avg <- getAverageLog(this, field="stdvs", indices=units, verbose=verbose);
  verbose && exit(verbose);

  medianSE <- getData(avg, indices=units, "intensities")$intensities;
  medianSE <- log2(medianSE);
  rm(avg);

  verbose && enter(verbose, "Calculating NUSE statistics for ", nbrOfArrays, 
                                                    " (specified) arrays");
  
  # For each file, calculate boxplot statistics
  stats <- list();
  for (kk in seq(along=arrays)) {
    array <- arrays[kk];
    cef <- getFile(this, array);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                          kk, getName(cef), nbrOfArrays));
    stdvs <- getData(cef, indices=units, "stdvs")$stdvs;
    stdvs <- log2(stdvs);
    stats[[kk]] <- boxplot.stats(stdvs/medianSE);
    rm(stdvs);
    verbose && exit(verbose);
  }
  rm(medianSE, units);
  names(stats) <- getNames(this)[arrays];
  verbose && exit(verbose);

  # Merge boxplot stats?
  if (merge) 
    stats <- mergeBoxplotStats(stats);

  attr(stats, "type") <- "NUSE";

  stats;
}) # calculateNuseStats()


##########################################################################
# HISTORY:
# 2008-02-25
# o Renamed to make it explicit that it is boxplot stats that are 
#   calculated.
# o Now calculate{Nuse|Rle}Stats() returns the list of boxplot stats,
#   not the combined ones.
# o Added boxplotStats() and support calculate{Nuse|Rle}BoxplotStats()
#   adopted from EPs code.
# 2008-02-22
# o Generalized from the QualityAssessmentModel.
##########################################################################
