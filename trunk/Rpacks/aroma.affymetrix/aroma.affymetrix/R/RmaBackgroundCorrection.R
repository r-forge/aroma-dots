###########################################################################/**
# @RdocClass RmaBackgroundCorrection
#
# @title "The RmaBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the RMA background adjustment function.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{addJitter}{}
#   \item{jitterSd}{}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
#*/###########################################################################
setConstructorS3("RmaBackgroundCorrection", function(..., addJitter=FALSE, jitterSd=0.2) {
  extend(BackgroundCorrection(..., typesToUpdate="pm"),
         "RmaBackgroundCorrection",
         .addJitter=addJitter,
         .jitterSd=jitterSd);
})


setMethodS3("getParameters", "RmaBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  if (this$.typesToUpdate=="pm") {
    pmonly <- TRUE;
  } else {
    pmonly <- FALSE;
  }
  
  # Get parameters of this class
  params2 <- list(
    addJitter = this$.addJitter,
    jitterSd = this$.jitterSd,
    pmonly = pmonly
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)

###########################################################################/**
# @RdocMethod process
#
# @title "Performs background correction"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "RmaBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Background correcting data set");

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected");
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

  # Get the output path
  outputPath <- getPath(this);
  mkdirs(outputPath);

  args <- c(list(ds, path=outputPath, verbose=verbose, overwrite=force), params);

  outputDataSet <- do.call("bgAdjustRma", args=args);
  
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Update the output data set
  this$.outputDataSet <- outputDataSet;

  verbose && exit(verbose);
  
  outputDataSet;
})



############################################################################
# HISTORY:
# 2007-03-21
# o Created.
############################################################################

