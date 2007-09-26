###########################################################################/**
# @RdocClass SrmaModel
#
# @title "The SrmaModel class"
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
setConstructorS3("SrmaModel", function(..., bandwidth=10e3, tags="*") {
  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(1,Inf));

  extend(ChromosomalModel(..., tags=tags), "SrmaModel",
    .outTuple = NULL,
    .shift = 0,
    .kernel = "gauss",
    .bandwidth = bandwidth
  )
})

setMethodS3("as.character", "SrmaModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", this, ...);
  s <- c(s, sprintf("Kernel: %s", this$.kernel));
  s <- c(s, sprintf("Bandwidth: %.2fkb", getBandwidth(this)/1e3));
  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)


setMethodS3("clearCache", "SrmaModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c()) {
    this[[ff]] <- NULL;
  }

  if (!is.null(this$.outTuple)) {
   this$.outTuple <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
})


setMethodS3("getAsteriskTag", "SrmaModel", function(this, ...) {
  classTag <- toupper(gsub("Model$", "", class(this)[1]));
  kernelTag <- this$.kernel;
  bandwidthTag <- sprintf("b=%d", getBandwidth(this));
  tags <- c(classTag, kernelTag, bandwidthTag);
#  tags <- paste(tags, collapse=",");
  tags;
}, protected=TRUE)


setMethodS3("getRootPath", "SrmaModel", function(this, ...) {
  tag <- getAsteriskTag(this)[1];
  sprintf("%sData", tolower(tag));
})

setMethodS3("getBandwidth", "SrmaModel", function(this, ...) {
  this$.bandwidth;
})

setMethodS3("setBandwidth", "SrmaModel", function(this, bandwidth, ...) {
  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(1,Inf));

  oldValue <- this$.bandwidth;
  this$.bandwidth <- bandwidth;
  invisible(oldValue);
})

setMethodS3("getOutputTuple", "SrmaModel", function(this, ..., force=FALSE,  vebose=FALSE) {
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
