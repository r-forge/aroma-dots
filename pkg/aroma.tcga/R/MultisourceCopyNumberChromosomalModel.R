###########################################################################/**
# @RdocClass MultisourceCopyNumberChromosomalModel
#
# @title "The MultisourceCopyNumberChromosomalModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a multi-source copy-number model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MultisourceChromosomalModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("MultisourceCopyNumberChromosomalModel", function(...) {
  extend(MultisourceChromosomalModel(...), "MultisourceCopyNumberChromosomalModel");
})




###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{See subclasses.}
# }
#
# \value{
#  See subclasses.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "MultisourceCopyNumberChromosomalModel", abstract=TRUE);




setMethodS3("extractRawCopyNumbers", "MultisourceCopyNumberChromosomalModel", function(this, ..., logBase=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getDouble(logBase, range=c(1, 10));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extract raw copy numbers");

  res <- extractRawGenomicSignals(this, ..., 
                  extractFcn=extractRawCopyNumbers, verbose=verbose);

  # Convert to the correct logarithmic base
  res <- extractRawCopyNumbers(res, logBase=logBase);

  verbose && print(verbose, res);
  verbose && exit(verbose);

  res;
}, protected=TRUE)



##############################################################################
# HISTORY:
# 2010-01-25
# o Created from CopyNumberChromosomalModel.
##############################################################################
