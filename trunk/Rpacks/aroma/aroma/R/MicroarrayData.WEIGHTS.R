#########################################################################/**
# @set "class=MicroarrayData"
# @RdocMethod setProbeWeights
#
# @title "Sets probe weights on one or several slides"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{slides}{The slides which should be assign weights.
#     If @NULL, all slides are included.}
#   \item{weights}{A @matrix or a @vector. If a @matrix it must have 
#     number of rows equal to number of spots, and number of columns equal 
#     to the number of specified slides.
#     If a @vector, it must be of length equal to the number of spots.
#     If @NULL, the weights are removed.}
# }
#
# \value{
#   Returns itself invisibly.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("setProbeWeights", "MicroarrayData", function(this, weights, slides=NULL, .force=FALSE) {
  if (!.force) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument: 'slides'
    slides <- validateArgumentSlides(this, slides=slides);
  
    # Argument 'weights'
    weights <- validateArgumentWeights(this, weights=weights, slides=slides);
    if (is.null(weights))
      throw("Argument 'weights' must not be NULL.");
  }

  if (!hasProbeWeights(this))
    resetProbeWeights(this);
  if (is.null(this$weights$probe)) {
    this$weights$probe <-  
                 matrix(1, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this));
  }

  if (is.null(slides)) {
    this$weights$probe <- weights;
  } else {
    this$weights$probe[,slides] <- weights;
  }
  
  invisible(this);
}) # setProbeWeights()


setMethodS3("getProbeWeights", "MicroarrayData", function(this, slides=NULL, missingValue=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);

  # Argument: 'missingValue'
  if (!is.null(missingValue)) {
    if (!is.numeric(missingValue) || length(missingValue) != 1 || 
        missingValue < 0 || missingValue > 1) {
      throw("Argument 'missingValue' must be a single value in [0,1].");
    }
  }

  res <- this$weights$probe[,slides];
  if (is.null(res)) {
    if (!is.null(missingValue))
      res <- matrix(missingValue, nrow=nbrOfSpots(this), ncol=length(slides));
  } else {
    res <- as.matrix(res);
  }

  res;
}) # getProbeWeights()


setMethodS3("resetProbeWeights", "MicroarrayData", function(this) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assure that all weight fields exist
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(this$weights))
    this$weights <- list();

  this$weights$probe <- NULL;

  invisible(this);
})


setMethodS3("hasProbeWeights", "MicroarrayData", function(this) {
  !is.null(this$weights$probe);
}) # hasProbeWeights()




setMethodS3("setSignalWeights", "MicroarrayData", function(this, weights, channel=NULL, slides=NULL, .force=FALSE) {
  if (!.force) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument: 'slides'
    slides <- validateArgumentSlides(this, slides=slides);
  
    # Argument: 'channels'
    channel <- validateArgumentChannel(this, channel=channel);
  
    # Argument: 'weights'
    weights <- validateArgumentWeights(this, weights=weights, slides=slides);
    if (is.null(weights))
      throw("Argument 'weights' must not be NULL.");
  }

  if (!hasSignalWeights(this, channel=channel))
    resetSignalWeights(this, channel=channel);

  if (is.null(this$weights$signal[[channel]])) {
    this$weights$signal[[channel]] <- 
                  matrix(1, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this));
  }
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(slides)) {
    this$weights$signal[[channel]] <- weights;
  } else {
    this$weights$signal[[channel]][,slides] <- weights;
  }
  
  invisible(this);
})



setMethodS3("getSignalWeights", "MicroarrayData", function(this, slides=NULL, channel=NULL, missingValue=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);

  # Argument: 'channel'
  channel <- validateArgumentChannel(this, channel=channel);

  # Argument: 'missingValue'
  if (!is.null(missingValue)) {
    if (!is.numeric(missingValue) || length(missingValue) != 1 || 
        missingValue < 0 || missingValue > 1) {
      throw("Argument 'missingValue' must be a single value in [0,1].");
    }
  }

  # Get the signal weights for the specified channel
  res <- this$weights$signal;
  if (!is.null(res)) {
    res <- res[[channel]];
    if (!is.null(res))
      res <- res[,slides, drop=FALSE];
  }

  if (is.null(res)) {
    if (!is.null(missingValue))
      res <- matrix(missingValue, nrow=nbrOfSpots(this), ncol=length(slides));
  } else {
    res <- as.matrix(res);
  }

  res;
}) # getSignalWeights()



setMethodS3("resetSignalWeights", "MicroarrayData", function(this, channel) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'channel'
  channel <- validateArgumentChannels(this, channel=channel);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assure that all weight fields exist
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(this$weights))
    this$weights <- list();

  if (is.null(this$weights$signal))
    this$weights$signal <- list();

  invisible(this);
})



setMethodS3("hasSignalWeights", "MicroarrayData", function(this, channel) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'channel'
  channel <- validateArgumentChannels(this, channel=channel);

  !is.null(this$weights$signal[[channel]]);
}) # hasSignalWeights()




setMethodS3("setWeights", "MicroarrayData", function(this, weights=NULL, .force=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.list(weights)) {
  } else if (!is.null(weights)) {
    throw("Unknown class of argument 'weights': ", class(weights)[1]);
  }

  this$weights <- weights;
}) # setWeights()



setMethodS3("getWeights", "MicroarrayData", function(this, slides=NULL) {
  weights <- list(
    probe = getProbeWeights(this, slides=slides)
  );
  weights;  
}) # getWeights()



setMethodS3("hasWeights", "MicroarrayData", function(this) {
  if (is.null(this$weights))
    return(FALSE);

  !all(sapply(this$weights, FUN=is.null))
}) # hasWeights()


setMethodS3("getWeightsAsString", "MicroarrayData", function(this) {
  w <- sapply(this$weights, FUN=is.null)
  if (all(w)) {
    "No weights are specified.";
  } else {
    paste("Weights (", paste(names(w), collapse=", "), 
          ") are specified.", sep="");
  }
}) # getWeightsAsString()



############################################################################
# HISTORY:
# 2005-02-14
# o BUG FIX: Used 'ch' instead of 'channel' in hasSignalWeights().
# 2005-02-09
# o Now all setNNNWeights() functions take argument '.force' to set the
#   weights without checking for consistency. Only for internal use!
# o Now all getNNNWeights() functions take argument 'missingValue' to
#   guarantee to get weights (default 1) even if missing.
# 2005-02-02
# o Now all weights types are stored under the field 'weights'.
# 2005-02-01
# o Added get- and setProbeWeights().
############################################################################
