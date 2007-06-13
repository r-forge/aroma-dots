###########################################################################/**
# @set "class=GladModel"
# @RdocMethod fitOne
#
# @title "Fits the GLAD model for one chromosome in one sample"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{ceList}{A @list of @see "ChipEffectFile" objects.}
#   \item{refList}{A @list of @see "ChipEffectFile" objects.}
#   \item{chromosome}{The chromosome for which the model should be fitted.}
#   \item{units}{(Optional) The subset of units to be matched.}
#   \item{useStdvs}{If @TRUE, standard deviations estimates for the 
#      chip effects are passed along and returned.  However, note that
#      GLAD is not making use of these estimates.}
#   \item{...}{Not used.}
#   \item{maxNAFraction}{A @numeric in [0,1] specifying the maximum fraction
#      of non-finite values allowed.  
#      If more are detected, this is interpreted as something has gone wrong
#      in the preprocessing and an error is thrown.}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @see "GLAD::profileCGH" object returned by @see "GLAD::glad".
# }
#
# @author
#
# \seealso{
#   Internally @see "GLAD::glad" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitOne", "GladModel", function(this, ceList, refList, chromosome, units=NULL, useStddvs=TRUE, ..., maxNAFraction=1/8, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Fitting GLAD");

  # Data set attributes
  chipTypes <- getChipTypes(this);
  arrayNames <- getArrays(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract arguments for glad().
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...);
  keep <- (names(args) %in% names(formals(GLAD::glad.profileCGH)));
  gladArgs <- args[keep];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x, M, stddev, chiptype, unit) from all chip types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving relative chip-effect estimates");
  # Get the chip types as a factor
  chipTypes <- as.factor(chipTypes);
  df <- NULL;
  for (kk in seq(along=chipTypes)) {
    chipType <- chipTypes[kk];
    verbose && enter(verbose, "Chip type: ", chipType);
    ce <- ceList[[kk]];
    if (!is.null(ce)) {
      ref <- refList[[kk]];
      df0 <- getXAM(ce, other=ref, chromosome=chromosome, units=units, verbose=less(verbose));
      df0 <- df0[,c("x", "M")];
      verbose && cat(verbose, "Number of units: ", nrow(df0));

      # Estimate the std dev of the raw log2(CN).  [only if ref is average across arrays]
      units0 <- as.integer(rownames(df0));
      # Get (mu, sigma) of theta (estimated across all arrays).
      data <- getDataFlat(ref, units=units0, verbose=less(verbose));
      # Number of arrays (for each unit)
      n <- readCel(getPathname(ref), indices=data$cell, readIntensities=FALSE, readPixels=TRUE)$pixels;
      # Use Gauss' approximation (since mu and sigma are on the intensity scale)
      sdM <- log2(exp(1)) * sqrt(1+1/n) * data$sdTheta / data$theta;
      rm(n);

      verbose && enter(verbose, "Scanning for non-finite values");
      n <- sum(!is.finite(df0[,"M"]));
      fraction <- n / nrow(df0);
      verbose && printf(verbose, "Number of non-finite values: %d (%.1f%%)\n", 
                                                             n, 100*fraction);
      if (fraction > maxNAFraction) {
        throw(sprintf("Something is wrong with the data. Too many non-finite values: %d (%.1f%% > %.1f%%)", as.integer(n), 100*fraction, 100*maxNAFraction));
      }
      verbose && exit(verbose);
  
      # Append SD, chip type, and CDF units.
      df0 <- cbind(df0, sdTheta=data$sdTheta, sdM=sdM, chipType=rep(chipType, length=length(units0)), unit=units0);
      rm(data);
  
      df <- rbind(df, df0);
      colnames(df) <- colnames(df0);
      rm(df0, units0);
    } else {
      verbose && cat(verbose, "No chip-effect estimates available sample: ", arrayNames[kk]);
    }

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    
    verbose && exit(verbose);
  } # for (kk in ...)
  

  verbose && enter(verbose, "Re-order by physical position");
  df <- df[order(df[,"x"]),];
  rownames(df) <- NULL;
  nbrOfUnits <- nrow(df);
  verbose && exit(verbose);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nbrOfUnits));
  verbose && exit(verbose);

  
  # Add T = M/sd (for future support to model (x,T) instead. /HB 2007-02-26)
  df <- cbind(df, T=df[,"M"]/df[,"sdM"]);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up GLAD data structure");
  # Put the data in a format recognized by GLAD
  df <- data.frame(
    LogRatio=df[,"M"], 
#   LogRatio=df[,"T"], 
    PosOrder=1:nbrOfUnits, 
    Chromosome=rep(chromosome, nbrOfUnits),
    PosBase=df[,"x"],
    # Add (chipType, units) identifiers to be able to backtrack SNP IDs etc.
    chipType=chipTypes[df[,"chipType"]],
    unit=df[,"unit"],
    # Add SD estimates
    sdTheta=df[,"sdTheta"],
    sdM=df[,"sdM"]
  );

  df <- GLAD::as.profileCGH(df);
  verbose && str(verbose, df);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling glad()");
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits);
  args <- c(list(df), gladArgs, list(verbose=as.logical(verbose)));
  rm(df, gladArgs);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # glad() writes to stdout; capture it and send it to the verbose object.
  stdout <- capture.output({
    fit <- do.call("glad", args);
  })
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;  
}, private=TRUE) # fitOne()


############################################################################
# HISTORY:
# 2007-06-12
# o Added argument 'maxNAFraction=1/8' to fitOne() of GladModel.
# 2007-06-11
# o Made explicit calls to GLAD::as.profileCGH() and GLAD::glad.profileCGH()
#   in fitOne() of GladModel.
# 2007-01-07
# o Renamed MultiGladModel to GladModel fully replacing the older class.
# o Added explicit garbage collection to minimize memory usage.
# o Added test against too (>10%) many non-finite values, which indicates
#   something is wrong with the theta estimates for either the sample or
#   the reference.
# 2006-12-21
# o BUG FIX: Made some updates to fitOne() so that the file-cache key would
#   be the same for all samples. Removed the caching completely, since fit()
#   now saves to file anyway.
# 2006-12-15
# o Works. Can now fit GLAD for multiple chip types together.
# o Created from deprecated fitGlad.CnChipEffectFile.R.
############################################################################
