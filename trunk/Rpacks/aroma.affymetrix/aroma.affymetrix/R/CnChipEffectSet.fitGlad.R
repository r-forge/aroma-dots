###########################################################################/**
# @set "class=CnChipEffectSet"
# @RdocMethod fitGlad
#
# @title "Fits the GLAD model to copy-number estimates"
#
# \description{
#  @get "title" of the files in the data set for a certain chromosome.
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
setMethodS3("fitGlad", "CnChipEffectSet", function(this, reference=NULL, arrays=1:nbrOfArrays(this), chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  throw("fitGlad() for CnChipEffectSet is deprecated since 2006-12-15.  Use the GladModel class instead.");

  require("GLAD") || throw("Package 'GLAD' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (is.null(reference)) {
  } else if (!inherits(reference, "CnChipEffectFile")) {
    throw("Argument 'reference' is not a CnChipEffectFile: ", 
                                                        class(reference)[1]);
  }

  # Argument 'chromosomes':
  chromosomes <- Arguments$getCharacters(chromosomes);
  chromosomes <- intersect(chromosomes, c(1:22, "X"));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting GLAD to data set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving reference chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(reference)) {
    verbose && enter(verbose, "No reference specified. Calculating average chip effects");
    reference <- getAverageFile(this, verbose=less(verbose));
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  for (chr in chromosomes) {
    res[[chr]] <- list();

    # Array by array
    for (aa in seq(along=arrays)) {
      array <- arrays[aa];
      ce <- getFile(this, array);
      verbose && enter(verbose, sprintf("Array %s (#%d of %d) on chromosome %s", 
                                      getName(ce), aa, length(arrays), chr));

      fit <- fitGlad(ce, reference=reference, chromosome=chr, ..., verbose=less(verbose));

      res[[chr]][[aa]] <- fit;

      rm(fit, ce, array);

      verbose && exit(verbose);
    } # for (aa in ...)

    verbose && enter(verbose, "Garbage collect");
    gc();
    verbose && exit(verbose);
  } # for (chr in ...)

  verbose && exit(verbose);

  res;  
}, private=TRUE, deprecated=TRUE) # fitGlad()


############################################################################
# HISTORY:
# 2007-06-11
# o BUG FIX: Called getFile(ces, ...) instead of  getFile(this, ...) in
#   fitGlad() of CnChipEffectSet.
# 2006-12-15
# o Made fitGlad() for CnChipEffectSet deprecated.  Use the GladModel class
#   instead.
# 2006-11-22
# o Updated fitGlad() so that a subset of chromosomes (and even arrays)
#   can be fitted.
# 2006-11-06
# o Created.
############################################################################
