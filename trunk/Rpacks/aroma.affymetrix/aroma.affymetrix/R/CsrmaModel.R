###########################################################################/**
# @RdocClass CsrmaModel
#
# @title "The CsrmaModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Chromosomal Smoothing Robust Multichip Analysis
#  method.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
#   \item{bandwidth}{A single @numeric specifying the smoothing bandwidth 
#     in units of nucleotides.}
#   \item{...}{Arguments passed to the constructor of 
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#*/###########################################################################
setConstructorS3("CsrmaModel", function(cesTuple=NULL, bandwidth=10e3, ...) {
  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(1,Inf));

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "CsrmaModel",
    .shift = 0,
    .bandwidth = bandwidth
  )
})


setMethodS3("getAsteriskTag", "CsrmaModel", function(this, ...) {
  "CSRMA";
}, protected=TRUE)

setMethodS3("getBandwidth", "CsrmaModel", function(this, ...) {
  this$.bandwidth;
})

setMethodS3("setBandwidth", "CsrmaModel", function(this, bandwidth, ...) {
  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(1,Inf));

  oldValue <- this$.bandwidth;
  this$.bandwidth <- bandwidth;
  invisible(oldValue);
})

setMethodS3("getOutputTuple", "CsrmaModel", function(this, ..., force=FALSE,  vebose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  outTuple <- this$.outTuple;
  if (force || is.null(outTuple)) {
    outTuple <- createOutputTuple(this, force=force, verbose=less(verbose, 2));
    this$.outTuple <- outTuple;
  }

  outTuple;
})


##############################################################################
# HISTORY:
# 2007-09-20
# o Created.
##############################################################################
