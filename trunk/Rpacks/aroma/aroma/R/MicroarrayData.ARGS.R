########################################################################/**
# @set "class=MicroarrayData"
# @RdocMethod validateArgumentSlide
#
# @title "Validates the argument 'slide'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{slide}{The argument value to be validated.}
# }
#
# \value{
#  Returns an interpreted and valid value. 
#  If the argument was invalid an @see "R.oo::Exception" is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("validateArgumentSlide", "MicroarrayData", function(this, slide=NULL) {
  arg <- slide;
  argStr <- paste(arg, collapse=", "); 

  # Interpret slide value
  if (is.null(slide)) {
    slide <- as.integer(1);
    if (nbrOfSlides(this) == 0)
      slide <- integer();
  } else if (identical(slide, "first")) {
    slide <- as.integer(1);
    if (nbrOfSlides(this) == 0)
      slide <- integer();
  } else if (identical(slide, "last")) {
    slide <- as.integer(nbrOfSlides(this));
    if (nbrOfSlides(this) == 0)
      slide <- integer();
  } else if (is.logical(slide)) {
    slide <- which(slide);
  } else if (is.character(slide)) {
    slide <- match(slide, getSlideNames(this), nomatch=-1);
  } else if (!is.numeric(slide)) {
    throw("Invalid value of argument 'slide': ", argStr);
  } 

  # Validate slide value
  if (length(slide) == 0) {
    throw("Argument 'slide' is empty or there is no such slide: ", argStr);
  } else if (length(slide) > 1) {
    throw("Argument 'slide' should only specify one single slide index: ", argStr);
  } else if (slide < 1 || slide > nbrOfSlides(this)) {
    throw("Argument 'slide' is out of range or non-existing: ", argStr);
  }

  as.integer(slide);
}, protected=TRUE)



#########################################################################/**
# @RdocMethod validateArgumentSlides
#
# @title "Validates the argument 'slides'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{The argument value to be validated.}
# }
#
# \value{
#  Returns an interpreted and valid value. 
#  If the argument was invalid an @see "R.oo::Exception" is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("validateArgumentSlides", "MicroarrayData", function(this, slides=NULL, unique=TRUE, sort=TRUE) {
  arg <- slides;
  argStr <- paste(arg, collapse=", "); 

  # Interpret slides value
  if (is.null(slides)) {
    slides <- seq(this);
    if (nbrOfSlides(this) == 0)
      return(integer());
  } else if (identical(slides, "all")) {
    slides <- seq(this);
  } else if (identical(slides, "first")) {
    slides <- as.integer(1);
    if (nbrOfSlides(this) == 0)
      slides <- integer();
  } else if (identical(slides, "last")) {
    slides <- as.integer(nbrOfSlides(this));
    if (nbrOfSlides(this) == 0)
      slides <- integer();
  } else if (identical(slides, "odd")) {
    slides <- seq(from=1, to=nbrOfSlides(this), by=2);
    if (nbrOfSlides(this) == 0)
      slides <- integer();
  } else if (identical(slides, "even")) {
    slides <- seq(from=2, to=nbrOfSlides(this), by=2);
    if (nbrOfSlides(this) < 2)
      slides <- integer();
  } else if (is.logical(slides)) {
    slides <- which(slides);
  } else if (is.numeric(slides)) {
    slidesIncl <- slides[slides >= 0];
    if (length(slidesIncl) == 0)
      slidesIncl <- seq(this);
    slidesExcl <- -slides[slides < 0];
    slides <- setdiff(slidesIncl, slidesExcl);
  } else if (is.character(slides)) {
    slides <- match(slides, getSlideNames(this), nomatch=-1);
  } else {
    throw("Invalid value of argument 'slides': ", argStr);
  }

  slides <- as.integer(slides);

  if (unique)
    slides <- unique(slides);
  if (sort)
    slides <- unique(slides);

  # Validate slide value
  if (length(slides) == 0) {
    throw("Argument 'slides' is empty or there are no such slides: ", argStr);
  } else if (any(slides < 1 | slides > nbrOfSlides(this))) {
    bad <- (slides < 1 | slides > nbrOfSlides(this));
    argStr <- paste(arg[bad], collapse=", ");
    throw("Argument 'slides' is out of range or non-existing: ", argStr);
  }

  as.integer(slides);
}, protected=TRUE)



setMethodS3("validateArgumentSpotIndex", "MicroarrayData", function(this, spotIndex=NULL, unique=TRUE, sort=TRUE) {
  arg <- spotIndex;
  argStr <- paste(arg, collapse=", "); 

  # Interpret spotIndex value
  if (is.null(spotIndex)) {
    spotIndex <- seq(length=nbrOfSpots(this));
  } else if (is.logical(spotIndex)) {
    spotIndex <- which(spotIndex);
  } else if (is.numeric(spotIndex)) {
    spotIndexIncl <- spotIndex[spotIndex >= 0];
    if (length(spotIndexIncl) == 0)
      spotIndexIncl <- seq(this);
    spotIndexExcl <- -spotIndex[spotIndex < 0];
    spotIndex <- setdiff(spotIndexIncl, spotIndexExcl);
  } else if (identical(spotIndex, "all")) {
    spotIndex <- seq(length=nbrOfSpots(this));
  } else if (identical(spotIndex, "first")) {
    spotIndex <- as.integer(1);
    if (nbrOfSpots(this) == 0)
      spotIndex <- integer();
  } else if (identical(spotIndex, "last")) {
    spotIndex <- as.integer(nbrOfSpots(this));
    if (nbrOfSpots(this) == 0)
      spotIndex <- integer();
  } else if (identical(spotIndex, "odd")) {
    spotIndex <- seq(from=1, to=nbrOfSpots(this), by=2);
    if (nbrOfSpots(this) == 0)
      spotIndex <- integer();
  } else if (identical(spotIndex, "even")) {
    spotIndex <- seq(from=2, to=nbrOfSpots(this), by=2);
    if (nbrOfSpots(this) < 2)
      spotIndex <- integer();
  } else {
    throw("Invalid value of argument 'spotIndex': ", argStr);
  }

  spotIndex <- as.integer(spotIndex);

  if (unique)
    spotIndex <- unique(spotIndex);
  if (sort)
    spotIndex <- unique(spotIndex);

  # Validate spotIndex value
  if (length(spotIndex) == 0) {
    throw("Argument 'spotIndex' is empty or there is no such spot: ", argStr);
  } else if (any(spotIndex < 1 | spotIndex > nbrOfSpots(this))) {
    bad <- (spotIndex < 1 | spotIndex > nbrOfSpots(this));
    argStr <- paste(arg[bad], collapse=", ");
    throw("Argument 'spotIndex' is out of range: ", argStr);
  }

  as.integer(spotIndex);
}, protected=TRUE)





#########################################################################/**
# @RdocMethod validateArgumentGroupBy
#
# @title "Validates the argument 'groupBy'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{groupBy}{Either a @see "LayoutGroups" object or a @character
#     string specifying a layout group. If @NULL, all spots are
#     considered.}
# }
#
# \value{
#  Returns a @list containing elements that each contains a @vector of
#  spot indices. If @NULL, the list will contain one element with 
#  \emph{all} indices.
#
#  If the argument was invalid an @see "R.oo::Exception" is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("validateArgumentGroupBy", "MicroarrayData", function(this, groupBy) {
  arg <- groupBy;
  argStr <- paste(arg, collapse=", "); 

  if (is.null(groupBy)) {
    spots <- list(1:nbrOfSpots(this));
  } else {
    # Identify the groups to be normalized seperately. 
    layout <- getLayout(this);
    if (is.null(layout))
      throw("Layout not specified.");
    if (!inherits(groupBy, "LayoutGroups")) {
      groupBy <- getLayoutGroupsByName(layout, groupBy);
    }
    spots <- getSpots(groupBy);
  }

  spots;
})


#########################################################################/**
# @RdocMethod validateArgumentWeights
#
# @title "Validates the argument 'weights'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{weights}{Either a @numeric, @logical @vector, or @NULL.
#     If a @logical, it is coersed to a zero-one @double @vector.
#     If a @numeric, each element must be in [0,1], otherwise an @Exception
#     is thrown.
#     If @NULL, then @NULL is returned.
#   }
#   \item{zeroOneOnly}{If @TRUE, only zero-one weights are allows. 
#     All non-zero weights are set to one.
#   }
# }
#
# \value{
#  Returns a @double @vector of weights of equal length as the input @vector.
#  If @NULL was given, then @NULL is returned.
#
#  If the argument was invalid an @see "R.oo::Exception" is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("validateArgumentWeights", "MicroarrayData", function(this, weights, slides=NULL, zeroOneOnly=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);

  arg <- weights;
  argStr <- paste(arg, collapse=", "); 

  if (is.null(weights)) {
    weights <- NULL;
  } else if (is.logical(weights)) {
    dim <- dim(weights);
    weights <- as.double(weights);
    dim(weights) <- dim;
  } else if (is.numeric(weights)) {
    dim <- dim(weights);
    weights <- as.double(weights);
    if (any(is.na(weights))) {
      throw("Argument 'weights' contains NA values.");
    }
    if (any(weights < 0 | weights > 1)) {
      throw("Argument 'weights' out of range [0,1]: ", 
           paste(weights[weights < 0.0 | weights > 1.0], collapse=", "));
    }
    if (zeroOneOnly && any(weights > 0 | weights < 1)) {
      weights <- as.double(as.integer(weights));
      warning("The current method only support zero-one {0,1} weights. Non-zero weights were set to one.");
    }
    dim(weights) <- dim;
  } else {
    throw("Argument 'weights' is of an unsupported datatype/class: ", class(weights)[1]);
  }

  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      if (nrow(weights) != nbrOfSpots(this))
        throw("The number of rows of argument 'weights' does not match the number of spots: ", nrow(weights));
      if (ncol(weights) != length(slides))
        throw("The number of columns of argument 'weights' does not match the number of slides: ", ncol(weights), " != ", nbrOfSlides(this));
    } else if (is.vector(weights)) {
      if (length(weights) != nbrOfSpots(this))
        throw("The length of argument 'weights' does not match the number of spots: ", length(weights));
      weights <- matrix(weights, nrow=nbrOfSpots(this), ncol=length(slides));
    }
  }

  weights;
})


setMethodS3("validateArgumentChannels", "MicroarrayData", function(this, channels=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Pre-process
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  knownChannels <- getChannelNames(this);
  nbrOfChannels <- length(knownChannels);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate channels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'channels'
  if (is.null(channels))
    return(knownChannels);

  if (is.numeric(channels)) {
    channels <- as.integer(channels);
    if (any(is.na(channels) | channels < 1 | channels > nbrOfChannels)) {
      throw("Argument 'channels' is out of range [1,", nbrOfChannels, "].");
    }
    channels <- knownChannels[channels];
  } else { 
    channels <- as.character(channels);
    channels <- toupper(channels);
    match <- match(channels, knownChannels);
    unknown <- is.na(match);
    if (any(unknown)) {
      throw("Argument 'channels' has unknown values: ", 
            paste(channels[unknown], collapse=", "));
    }
  }
  channels;
}) # validateArgumentChannels()


setMethodS3("validateArgumentChannel", "MicroarrayData", function(this, channel=NULL) {
  channel <- validateArgumentChannels(this, channels=channel);
  if (length(channel) != 1)
    throw("Argument 'channel' must be specify exactly one channel, not ", 
                                                           length(channel));
  channel;
})



############################################################################
# HISTORY:
# 2005-02-28
# o BUG FIX: validateArgumentGroupBy() was looking for 'LayoutGroup' and 
#   not 'LayoutGroups'.
# 2005-02-09
# o Now validateArgumentChannels() uses getChannelNames() and is therefore
#   general for all MicroarrayData classes overriding getChannelNames().
# 2005-02-08
# o Added validateArgumentChannels() for RGData and 
#   validateArgumentChannel() for MicroarrayData.
# 2005-02-01
# o Replaced all stop() with throws().
# 2005-01-23
# o Some simple tests have been run on validateArgumentWeights().
# 2004-12-27
# o Added validateArgumentWeights().
# 2004-12-12
# o Added validateArgumentGroupBy().
# 2004-12-01
# o BUG FIX: Using slide names with  validateArgumentSlides() gave an
#   error due to a typo.
# 2003-12-30
# o Added support for slide names in validateArgumentSlide[s]().
# 2003-10-30
# o Created. All methods should validate and parse the arguments in the
#   same way. This will result in a more consistent user-interface, but also
#   more robust design and code. 
#   "A small step for mankind a big step for man..."
############################################################################
