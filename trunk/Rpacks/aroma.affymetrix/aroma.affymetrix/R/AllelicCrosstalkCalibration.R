###########################################################################/**
# @RdocClass AllelicCrosstalkCalibration
#
# @title "The AllelicCrosstalkCalibration class"
#
# \description{
#  @classhierarchy
#
#  This class represents a calibration function that transforms the 
#  probe-level signals such that the signals from the two alleles are 
#  orthogonal.
#  The method fits and calibrates PM signals only.  MM signals will not 
#  affect the model fitting and are unaffected.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#     @see "ProbeLevelTransform".}
#   \item{targetAvg}{The signal that average allele A and average allele B
#     signals should have after calibration.}
#   \item{alpha, q, Q}{}
# }
#
# \section{What probe signals are updated}{ 
#    Calibration for crosstalk between allele signals applies by definition
#    only SNP units.
#    Note that, non-calibrated signals will be saved in the output files.
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AllelicCrosstalkCalibration", function(..., targetAvg=2200, alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98) {
  extend(ProbeLevelTransform(...), "AllelicCrosstalkCalibration",
    .targetAvg = targetAvg,
    .alpha = alpha,
    .q = q,
    .Q = Q
  )
})



setMethodS3("getParameters", "AllelicCrosstalkCalibration", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  defaults <- formals(calibrateAllelicCrosstalk.AffymetrixCelFile);
  params$targetAvg <- this$.targetAvg;
  params$alpha <- this$.alpha;
  params$q <- this$.q;
  params$Q <- this$.Q;

  params;
}, private=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Calibrates the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already calibrated is re-calibrated, 
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @examples "../incl/normalizeQuantile.Rex"
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "AllelicCrosstalkCalibration", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calibrates data set for allelic cross talk");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already calibrated");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
  outputPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- c(list(ds, path=outputPath, verbose=verbose), params);
  rm(params);

  outputDataSet <- do.call("calibrateAllelicCrosstalk", args=args);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Update the output data set
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);
  
  outputDataSet;
})

############################################################################
# HISTORY:
# 2007-03-29
# o Now 'targetAvg' defaults to 2200 so that allele A and allele B signals
#   are rescaled to be one the same scale.  If so,  it does not make sense
#   to do background correction afterwards.
# o Added getParameters().
# o Added support for arguments 'targetAvg', 'alpha', 'q', and 'Q'.
# 2006-12-08
# o Now this class inherits from the ProbePreprocessing class.
# o Now this pre-processor output results to probeData/.
# o Renamed from AllelicCrosstalkCalibrator.
# 2006-11-18
# o Removed version and subversion tags, and related functions. 
#   Now getTags() returns the tags of the input data set plus any tags 
#   of this instance.
# 2006-11-02
# o Created from QuantileNormalizer.R.
############################################################################
