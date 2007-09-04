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
#   \item{data}{A @data.frame with columns \code{M} (log-ratio) and 
#      \code{x} (locus position).
#   }
#   \item{chromosome}{An @integer specifying the index of the chromosome to
#      be fitted.}
#   \item{...}{Additional arguments passed down to the internal fit function.}
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
setMethodS3("fitOne", "GladModel", function(this, data, chromosome, ..., verbose=FALSE) {
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract arguments for glad().
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...);
  keep <- (names(args) %in% names(formals(GLAD::glad.profileCGH)));
  fitArgs <- args[keep];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up GLAD data structure");
  # Put the data in a format recognized by GLAD
  nbrOfUnits <- nrow(data);
  chipTypes <- getChipTypes(this);
  data <- data.frame(
    LogRatio=data[,"M"], 
#   # Use T = M/sd (for modelling (x,T) instead. /HB 2007-02-26)
#   LogRatio=data[,"M"]/data[,"sdM"], 
    PosOrder=1:nbrOfUnits, 
    Chromosome=rep(chromosome, nbrOfUnits),
    PosBase=data[,"x"],
    # Add (chipType, units) identifiers to be able to backtrack SNP IDs etc.
    chipType=chipTypes[data[,"chipType"]],
    unit=data[,"unit"],
    # Add SD estimates
    sdTheta=data[,"sdTheta"],
    sdM=data[,"sdM"]
  );
  data <- GLAD::as.profileCGH(data);
  verbose && str(verbose, data);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling glad()");
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits);
  args <- c(list(data), fitArgs, list(verbose=as.logical(verbose)));
  rm(data, fitArgs);

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
# 2007-08-20
# o Initial tests show that the updated GladModel gives identical results.
# o Moved the code to extract raw CN data to superclass, see getRawCnData().
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
