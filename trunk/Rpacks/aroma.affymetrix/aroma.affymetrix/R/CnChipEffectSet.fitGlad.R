###########################################################################/**
# @set "class=CnChipEffectSet"
# @RdocMethod fitGlad
#
# @title "Fits the GLAD model to copy-number estimates"
#
# \description{
#  @get "title" of all files in the data set for a certain chromosome.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to \code{fitGlad()} of 
#     @see "CnChipEffectFile".}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns...
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitGlad", "CnChipEffectSet", function(this, reference=NULL, ..., verbose=FALSE) {
  require(GLAD) || throw("Package 'GLAD' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (is.null(reference)) {
  } else if (!inherits(reference, "CnChipEffectFile")) {
    throw("Argument 'reference' is not a CnChipEffectFile: ", 
                                                        class(reference)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting GLAD to data set");

  if (is.null(reference)) {
    reference <- getAverageFile(this, verbose=less(verbose));
  }

  fit <- lapply(this, fitGlad, reference=reference, ..., verbose=less(verbose));
  verbose && exit(verbose);

  fit;  
}) # fitGlad()


############################################################################
# HISTORY:
# 2006-11-06
# o Created.
############################################################################
