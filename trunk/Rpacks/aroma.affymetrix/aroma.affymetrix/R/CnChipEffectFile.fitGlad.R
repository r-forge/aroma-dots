###########################################################################/**
# @set "class=CnChipEffectFile"
# @RdocMethod fitGlad
#
# @title "Fits the GLAD model to copy-number estimates"
#
# \description{
#  @get "title" of a certain chromosome.
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The @see "CnChipEffectFile" object used as the 
#     reference chip effects when calculating the raw (relative) copy-number
#     estimates..}
#   \item{chromosome}{The chromosome for which the model should be fitted.}
#   \item{units}{(Optional) The subset of units to be matched.}
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
setMethodS3("fitGlad", "CnChipEffectFile", function(this, reference, chromosome, units=NULL, ..., force=FALSE, verbose=FALSE) {
  require(GLAD) || throw("Package 'GLAD' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (!inherits(reference, "CnChipEffectFile")) {
    throw("Argument 'reference' is not a CnChipEffectFile: ", 
                                                        class(reference)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting GLAD");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="fitGlad", sample=getIdentifier(this), reference=getIdentifier(reference), chromosome=chromosome, units=units, ...);
  fit <- loadCache(key=key);
  if (!is.null(fit) && !force) {
    verbose && cat(verbose, "Cached in memory.");
    verbose && exit(verbose);
    return(fit);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (X, M, A)
  verbose && enter(verbose, "Retrieving relative copy-number estimates");
  df <- getXAM(this, other=reference, chromosome=chromosome, units=units, verbose=less(verbose));
  verbose && enter(verbose, "Re-order by physical position");
  df <- df[order(df[,"x"]),];
  verbose && exit(verbose);
  verbose && str(verbose, df);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nrow(df)));

  # Put the data in a format recognized by GLAD
  df <- data.frame(
    LogRatio=unname(df[,"M"]), 
    PosOrder=1:nrow(df), 
    Chromosome=rep(chromosome, nrow(df)), 
    PosBase=unname(df[,"x"])
  );
  verbose && str(verbose, df);
  df <- as.profileCGH(df);
  verbose && str(verbose, df);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling glad()");
  fit <- glad(df, ..., verbose=as.logical(verbose));
  verbose && exit(verbose);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  comment <- paste(unlist(key), collapse=";");
  saveCache(fit, key=key, comment=comment);

  fit;  
}) # fitGlad()


############################################################################
# HISTORY:
# 2006-10-30
# o Added Rdoc comments.
# o BUG FIX: glad() requires data points to be order by physical position,
#   which was not (necessarily) the case.
# o Added more argument validation.
# 2006-10-17
# o Created.
############################################################################
