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
#   \item{ceList}{A @list of @see "CnChipEffectFile" objects.}
#   \item{refList}{A @list of @see "CnChipEffectFile" objects.}
#   \item{chromosome}{The chromosome for which the model should be fitted.}
#   \item{units}{(Optional) The subset of units to be matched.}
#   \item{useStdvs}{If @TRUE, standard deviations estimates for the 
#      chip effects are passed along and returned.  However, note that
#      GLAD is not making use of these estimates.}
#   \item{...}{Not used.}
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
setMethodS3("fitOne", "MultiGladModel", function(this, ceList, refList, chromosome, units=NULL, useStddvs=TRUE, ..., force=FALSE, verbose=FALSE) {
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


  fullname <- getFullName(this);
  chipTypes <- getChipTypes(this);
  arrayNames <- getArrays(this);

  ceNames <- sapply(ceList, FUN=getFullName);
  refNames <- sapply(refList, FUN=getFullName);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract arguments for glad().
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...);
  keep <- (names(args) %in% names(formals(glad.profileCGH)));
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
  
      verbose && enter(verbose, "Retrieving stddvs of chip effects");
      units0 <- as.integer(rownames(df0));
      sdTheta <- getDataFlat(ce, units=units0, fields="sdTheta", verbose=less(verbose))[,"sdTheta"];
      verbose && exit(verbose);
  
      # Append SD, chip type, and CDF units.
      df0 <- cbind(df0, sdTheta, chipType=rep(chipType, length=length(units0)), unit=units0);
      rm(sdTheta);
  
      df <- rbind(df, df0);
      colnames(df) <- colnames(df0);
      rm(df0, units0);
    } else {
      verbose && cat(verbose, "No chip-effect estimates available sample: ", arrayNames[kk]);
    }
    verbose && exit(verbose);
  }
  
  verbose && enter(verbose, "Re-order by physical position");
  df <- df[order(df[,"x"]),];
  rownames(df) <- NULL;
  nbrOfUnits <- nrow(df);
  verbose && exit(verbose);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nbrOfUnits));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Put the data in a format recognized by GLAD
  df <- data.frame(
    LogRatio=df[,"M"], 
    PosOrder=1:nbrOfUnits, 
    Chromosome=rep(chromosome, nbrOfUnits),
    PosBase=df[,"x"],
    # Add (chipType, units) identifiers to be able to backtrack SNP IDs etc.
    chipType=chipTypes[df[,"chipType"]],
    unit=df[,"unit"],
    # Add SD estimates
    sdTheta=df[,"sdTheta"]
  );
  df <- as.profileCGH(df);
  verbose && str(verbose, df);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling glad()");
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits);
  args <- c(list(df), gladArgs, list(verbose=as.logical(verbose)));
  fit <- do.call("glad", args);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;  
}) # fitOne()


############################################################################
# HISTORY:
# 2006-12-21
# o BUG FIX: Made some updates to fitOne() so that the file-cache key would
#   be the same for all samples. Removed the caching completely, since fit()
#   now saves to file anyway.
# 2006-12-15
# o Works. Can now fit GLAD for multiple chip types together.
# o Created from deprecated fitGlad.CnChipEffectFile.R.
############################################################################
