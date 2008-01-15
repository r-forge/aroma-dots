#########################################################################/**
# @RdocClass MicroarrayData
#
# @title "The MicroarrayData class"
#
# \description{
#  @classhierarchy
#
#  This class is \emph{abstract} can not be instanciated.
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#    \tab \code{layout} \tab The layout structure of the microarray slide(s). \cr
#  }
#
#  @allmethods
# }
#
# @author
#*/#########################################################################
setConstructorS3("MicroarrayData", function(layout=NULL, extras=list()) {
  if (!is.null(layout)) {
    if (!inherits(layout, "Layout"))
      throw("Argument 'layout' is not of class Layout: ", data.class(layout));
  }
  
  extend(Object(), "MicroarrayData", 
    layout       = layout,
    .isSaturated = list(),
    .treatments  = NULL,
    .fieldNames  = NULL,
    .exclude     = list(NULL),
    .extras      = extras,
    .labels      = list(),
    .flags       = list()
  )
}, abstract=TRUE);



############################################################################
############################################################################
## 
##  DISPLAYING (textual) AND SUMMARIZING METHODS
## 
############################################################################
############################################################################


setMethodS3("as.character", "MicroarrayData", function(this) {
  if (isAbstract(Class$forName(class(this)[1])))
    throw("Class is abstract: ", class(this)[1]);

  s <- paste(sep="", class(this)[1], ":");
  s <- paste(sep="", s, " Number of slides: ", nbrOfSlides(this), ".");
  fields <- getFieldNames(this);
  if (length(fields) <= 5) {
    fields <- paste(fields, collapse=", ");
    s <- paste(s, " Fields: ", fields, ".", sep="")
  } else
    s <- paste(s, " Number of fields: ", length(fields), ".", sep="");
  if (hasLayout(this))
    s <- paste(sep="", s, " ", getLayout(this));
  if (hasWeights(this)) {
    w <- sapply(this$weights, FUN=is.null);
    s <- paste(s, " Weights (", paste(names(w), collapse=", "), 
                  ") are specified.", sep="");
  }
  if (!is.null(this$.cache) && length(this$.cache) > 0)
    s <- paste(sep="",s,", cached values for ",
               paste(names(this$.cache), collapse=", "), " exist");
  s;
})







#########################################################################/**
# @RdocMethod str
#
# @title "Compactly Display the Structure of a MicroarrayData object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# @author
#
# \seealso{
#   @see "base::str".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("str", "MicroarrayData", function(object) {
  # To please R CMD check...
  this <- object;

  str(getLayout(this));
  str(this$.repfcn);
  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
  for (field in getFieldNames(this))
    str(this[[field]]);
})





#########################################################################/**
# @RdocMethod summary
#
# @title "Gets summary statistics for the fields"
#
# \description{
#   @get "title". 
#   This is identical to calling @see "base::summary" on the result from
#   @seemethod "as.data.frame".
# }
#
# @synopsis
#
# @author
#
# \seealso{
#   @seemethod "as.data.frame".
#   @see "base::summary".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("summary", "MicroarrayData", function(object, ...) {
  # To please R CMD check...
  this <- object;

  summary(as.data.frame(this), ...);
})




############################################################################
############################################################################
## 
##  ADDING, REMOVING AND KEEPING SLIDES
## 
############################################################################
############################################################################


#########################################################################/**
# @RdocMethod getChannelNames
#
# @title "Gets the names of the channels"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{digitsOnly}{If @TRUE, all non-digit characters are excluded from
#     the channel names.}
#   \item{settings}{Internal use only.}
# }
#
# \value{
#  Returns an @character string @vector. An element with value @NA is a
#  channel that had a zero-length name.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################s
setMethodS3("getChannelNames", "MicroarrayData", function(this) {
  return(NULL);
})



#########################################################################/**
# @RdocMethod append
#
# @title "Appends another MicroarrayData object"
#
# @synopsis
#
# \arguments{
#  \item{...}{Other MicroarrayData object(s) to be appended to this object.}
# }
#
# \description{
#   Appends another identically structured MicroarrayData object(s).
# }
#
# \examples{\dontrun{
#   # Example that reads all GenePixData files in the current directory
#   # and extracts the log-ratios and the log-intensities into *one*
#   # MAData object. This approach only stores one GPR file in memory
#   # at the time, which makes it possible to read in approximately
#   # 10-20 times more files.
#   files <- list.files(pattern="*.gpr");
#   ma <- MAData();
#   for (file in files) {
#     gpr <- GenePixData$read(file);
#     raw <- getRawData(gpr);
#     rm(gpr); # Delete as quick as possible to optimize memory.
#     append(ma, getSignal(raw, bgSubtract=TRUE));
#     rm(raw); # Delete as quick as possible to optimize memory.
#   }
# }}
#
# @author
#
# \seealso{
#   @seemethod "removeSlides" and @seemethod "keepSlides".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("append", "MicroarrayData", function(this, other) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!inherits(other, "MicroarrayData"))
    throw("Argument 'other' is not of class MicroarrayData: ", data.class(other));
  if (data.class(this) != data.class(other))
    warning("Trying to append a MicroarrayData object of a different class.");

  # Assert that the two MicroarrayData objects have the same array layout.
  tLayout <- getLayout(this);
  oLayout <- getLayout(other);
  if (inherits(oLayout, "Layout")) {
    if (is.null(tLayout)) {
      # i) If 'this' does *not* have a layout set, set it...
      setLayout(this, oLayout);
    } else if (!equals(oLayout, tLayout)) {
      # ii) ...otherwise, assert they have the same layout.
      throw("Can not append MicroarrayData object with a different layout.");
    }
  }

  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
  attr(other, ".memberAccessorOrder") <- c(2,3,1,4,5);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Probe weights
  if (hasProbeWeights(this) || hasProbeWeights(other)) {
    w <- cbind(
      getProbeWeights(this, missingValue=1),
      getProbeWeights(other, missingValue=1)
    );
    setProbeWeights(this, weights=w, .force=TRUE);
    rm(w);
  }

  # Signal weights
  for (ch in getChannelNames(this)) {
    if (hasSignalWeights(this, channel=ch) || 
        hasSignalWeights(other, channel=ch)) {
      w <- cbind(
        getSignalWeights(this, channel=ch, missingValue=1),
        getSignalWeights(other, channel=ch, missingValue=1)
      );
      setSignalWeights(this, channel=ch, weights=w, .force=TRUE);
      rm(w);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Other fields
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fields <- intersect(getFieldNames(this), getFieldNames(other));
  for (field in fields) {
    this[[field]] <- cbind(this[[field]], other[[field]]);
    otherSaturated <- other$.isSaturated[[field]];
    if (is.null(otherSaturated))
      otherSaturated <- matrix(NA, nrow=nbrOfSpots(this), ncol=ncol(other[[field]]));
    this$.isSaturated[[field]] <- cbind(this$.isSaturated[[field]], otherSaturated);
  }

  invisible(this);
})



#########################################################################/**
# @RdocMethod removeSlides
#
# @title "Remove specified slides"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{slides}{@vector of slides to be removed. 
#    If @NULL, all slides are removed.}
# }
#
# @author
#
# \seealso{
#   @seemethod "append" and @seemethod "keepSlides".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("removeSlides", "MicroarrayData", function(this, slides) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);

  for (field in getFieldNames(this)) {
    this[[field]] <- this[[field]][,-slides, drop=FALSE];
    this$.isSaturated[[field]] <- this$.isSaturated[[field]][,-slides, drop=FALSE];
  }

  # Weights
  if (length(this$weights) > 0) {
    # Signal weights
    if (length(this$weights$signal) > 0) {
      for (ch in 1:length(this$weights$signal)) {
        this$weights$signal[[ch]] <- 
                             this$weights$signal[[ch]][,-slides, drop=FALSE];
      }
    }

    # Probe weights
    if (length(this$weights$probe) > 0) {
      this$weights$probe <- this$weights$probe[,-slides, drop=FALSE];
    }
  }
})




#########################################################################/**
# @RdocMethod keepSlides
#
# @title "Remove all but the specified slides"
#
# \description{
#   @get "title".
#
#   Note that this method can also be used to reshuffle the slides in a
#   certain order (and even duplicate some slides).
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{@vector of slides to be kept. If @NULL, nothing is done.}
# }
#
# @author
#
# \seealso{
#   @seemethod "append" and @seemethod "removeSlides".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("keepSlides", "MicroarrayData", function(this, slides) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);

  for (field in getFieldNames(this)) {
    this[[field]] <- this[[field]][,slides, drop=FALSE];
    this$.isSaturated[[field]] <- 
                            this$.isSaturated[[field]][,slides, drop=FALSE];
  }

  # Weights
  if (length(this$weights) > 0) {
    # Signal weights
    if (length(this$weights$signal) > 0) {
      for (ch in 1:length(this$weights$signal)) {
        this$weights$signal[[ch]] <- 
                             this$weights$signal[[ch]][,slides, drop=FALSE];
      }
    }

    # Probe weights
    if (length(this$weights$probe) > 0) {
      this$weights$probe <- this$weights$probe[,slides, drop=FALSE];
    }
  }
})



#########################################################################/**
# @RdocMethod keepSpots
#
# @title "Remove all but the specified spots from all slides"
#
# \description{
#   @get "title".
#
#   Note that this method can also be used to reshuffle the spots in a
#   certain order (and even duplicate some spots).
# }
#
# @synopsis
#
# \arguments{
#   \item{index}{@vector of spot indices to be kept. If @NULL, nothing
#    is done.}
# }
#
# @author
#
# \seealso{
#   @seemethod "removeSpots".
#   For similar operations on slides, see @seemethod "keepSlides" and
#   @seemethod "removeSlides".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("keepSpots", "MicroarrayData", function(this, index) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'index'
  if (any(index < 1 | index > nbrOfSpots(this)))
    throw("Argument 'index' out of range.");

  for (field in getFieldNames(this)) {
    this[[field]] <- this[[field]][index,];
  }

  # Weights
  if (length(this$weights) > 0) {
    # Signal weights
    if (length(this$weights$signal) > 0) {
      for (ch in 1:length(this$weights$signal)) {
        this$weights$signal[[ch]] <- 
                              this$weights$signal[[ch]][index,, drop=FALSE];
      }
    }

    # Probe weights
    if (length(this$weights$probe) > 0) {
      this$weights$probe <- this$weights$probe[index,, drop=FALSE];
    }
  }
})




#########################################################################/**
# @RdocMethod removeSpots
#
# @title "Remove specified spots"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{spots}{@vector of spots to be removed.
#               If @NULL, all spots are removed.}
# }
#
# @author
#
# \seealso{
#   @seemethod "keepSpots".
#   For similar operations on slides, see @seemethod "keepSlides" and
#   @seemethod "removeSlides".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("removeSpots", "MicroarrayData", function(this, index) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'index'
  if (any(index < 1 | index > nbrOfSpots(this)))
    throw("Argument 'index' out of range.");

  for (field in getFieldNames(this)) {
    this[[field]] <- this[[field]][-index,];
  }

  # Weights
  if (length(this$weights) > 0) {
    # Signal weights
    if (length(this$weights$signal) > 0) {
      for (ch in 1:length(this$weights$signal)) {
        this$weights$signal[[ch]] <- 
                              this$weights$signal[[ch]][-index,, drop=FALSE];
      }
    }

    # Probe weights
    if (length(this$weights$probe) > 0) {
      this$weights$probe <- this$weights$probe[-index,, drop=FALSE];
    }
  }
})




#########################################################################/**
# @RdocMethod isSaturated
#
# @title "Check if some observations are saturated or not"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{channels}{@vector of channels to be checked. If @NULL, all are checked.}
#   \item{index}{@vector of spot indicies to be checked. If @NULL, all are checked.}
#   \item{slides}{@vector of slides to be checked. If @NULL, all are checked.}
#   \item{na.rm}{If @TRUE, observations for which it is unknown (==@NA) if 
#     it is saturated or not, are excluded.}
# }
#
# @author
#
# \seealso{
#   @seemethod "setSaturated".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("isSaturated", "MicroarrayData", function(this, channels=NULL, index=NULL, slides=NULL, na.rm=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(index)) {
    index <- 1:nbrOfSpots(this);
  } else if (is.logical(index)) {
    if (length(index) > nbrOfSpots(this))
      throw("Argument 'index' with logical values is too long: ", length(index));
    index <- which(index);
  } else if (is.numeric(index)) {
    if (any(index < 1 || index > nbrOfSpots(this)))
      throw("Argument 'index' is out of range.");
  } else {
    throw("Unknown value(s) of argument 'index'.");
  }

  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(channels)) {
    channels <- getFieldNames(this);
  } else if (!all(channels %in% getFieldNames(this))) {
    throw("Unknown values in argument 'channels': ", paste(channels, collapse=", "));
  }

  result <- NULL;  
  for (channel in channels) {
    value <- this$.isSaturated[[channel]][index, slides];
    if (!is.null(value)) {
      if (is.null(result)) {
        result <- value;
      } else {
        if (na.rm == TRUE)
          value[is.na(value)] <- FALSE;
        result <- result |  value;
      }
    }
  }

  result;
}, private=TRUE, trial=TRUE)


#########################################################################/**
# @RdocMethod setSaturated
#
# @title "Set which observations in a given channel that are saturated or not"
#
# \description{
#   @get "title" (or unknown).
# }
#
# @synopsis
#
# \arguments{
#   \item{channels}{@vector of channels to be set. 
#     If @NULL, all are considered.}
#   \item{index}{@vector of spot indicies to be set. 
#     If @NULL, all are considered.}
#   \item{slides}{@vector of slides to be set.
#     If @NULL, all are considered.}
#   \item{status}{A @matrix of size matching the \code{index} and 
#     \code{slides} arguments (if not the values will be looped over).
#     If a value is @TRUE, the corresponding observation is saturated. 
#     If @FALSE, it is not and if @NA, it is unknown. }
# }
#
# @author
#
# \seealso{
#   @seemethod "isSaturated".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setSaturated", "MicroarrayData", function(this, channel, index=NULL, slides=NULL, status=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.element(channel, getFieldNames(this)))
    throw("Unknown channel: ", channel);

  if (is.null(index)) {
    if (length(index) > nbrOfSpots(this))
      throw("Argument 'index' with logical values is too long: ", length(index));
    index <- 1:nbrOfSpots(this);
  } else if (is.numeric(index)) {
    if (any(index < 1 || index > nbrOfSpots(this)))
      throw("Argument 'index' is out of range.");
  } else {
    throw("Unknown value(s) of argument 'index'.");
  }

  slides <- validateArgumentSlides(this, slides=slides);
  
  oldStatus <- this$.isSaturated[[channel]];

  status <- matrix(status, nrow=length(index), ncol=length(slides));
  if (is.null(oldStatus)) {
    this$.isSaturated[[channel]] <- matrix(FALSE, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this));
  }
  this$.isSaturated[[channel]][index,slides] <- status;

  invisible(oldStatus);
}, private=TRUE, trial=TRUE)


############################################################################
############################################################################
## 
##  TEST AND COMPARISON METHODS
## 
############################################################################
############################################################################



###########################################################################/**
# @RdocMethod equals
#
# @title "Compares a MicroarrayData object with another Object"
#
# @synopsis
#
# \description{
#  @get "title" and returns
#  @TRUE if they are equal. For the MicroarrayData object \code{raw}
#  to be equal to another Object \code{obj} it must be true that:
#
#  \item{1}{\code{obj} is of class MAData,}
#  \item{2}{\code{size(raw) == size(obj)},}
#  \item{3}{\code{equals(getLayout(raw), getLayout(obj))} is true, and}
#  \item{4}{the Euclidian distance between the each field in the two objects
#    is small.}
# }
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   ma2 <- clone(ma)
#   print(equals(ma, ma2))   # TRUE
#
#   layout2 <- Layout(4,4,21,19)
#   setLayout(ma2, layout2)
#   print(equals(ma, ma2))   # FALSE
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("equals", "MicroarrayData", function(this, other) {
  res <- FALSE;

  if (identical(this, other))
    return(TRUE);
  
  if (!inherits(other, "MicroarrayData")) {
    attr(res, "reason") <- 
     paste("Other object is not a MicroarrayData object: ", class(other)[1]);
    return(res);
  }

  if (size(this) != size(other)) {
    attr(res, "reason") <- "Sizes differ";
    return(res);
  }

  if (!equals(getLayout(this), getLayout(other))) {
    attr(res, "reason") <- "Array layouts differ";
    return(res);
  }

  fieldNames <- getFieldNames(this);
  if (!all(is.element(fieldNames, getFieldNames(other)))) {
    attr(res, "reason") <- "Objects have different sets of fields.";
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Field names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  eps <- 1e-12*size(this);
  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
  attr(other, ".memberAccessorOrder") <- c(2,3,1,4,5);
  for (field in fieldNames) {
    if (is.numeric(this[[field]])) {
      if ( sum((this[[field]]-other[[field]])^2, na.rm=TRUE) > eps ) {
        attr(res, "reason") <- 
                       paste("numeric field '", field, "' differ", sep="");
        return(res);
      }
    } else if (this[[field]] != other[[field]]) {
      res <- FALSE;
      attr(res, "reason") <- paste("field '", field, "' differ", sep="");
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Probe weights
  if (hasProbeWeights(this) || hasProbeWeights(other)) {
    w1 <- getProbeWeights(this, missingValue=1);
    w2 <- getProbeWeights(other, missingValue=1);
    dw <- sum((w1-w2)^2, na.rm=TRUE);
    rm(w1,w2);
    if (dw > eps) {
      res <- FALSE;
      attr(res, "reason") <- "probe weights differ";
      return(res);
    }
  }

  # Signal weights
  for (ch in getChannelNames(this)) {
    if (hasSignalWeights(this, channel=ch) || 
        hasSignalWeights(other, channel=ch)) {
      w1 <- getSignalWeights(this, channel=ch, missingValue=1);
      w2 <- getSignalWeights(other, channel=ch, missingValue=1);
      dw <- sum((w1-w2)^2, na.rm=TRUE);
      rm(w1,w2);
      if (dw > eps) {
        res <- FALSE;
        attr(res, "reason") <- 
                  paste("signal weights in channel '", ch, "' differ", sep="");
        return(res);
      }
    }
  }
  
  TRUE;
})




setMethodS3("getExcludedSpots", "MicroarrayData", function(this, slides=NULL) {
  if (!hasExcludedSpots(this))
    matrix(FALSE, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this))
  else
    unwhich(this$.exclude, dim=c(nbrOfSpots(this), nbrOfSlides(this)));
}, trial=TRUE)


setMethodS3("hasExcludedSpots", "MicroarrayData", function(this, slides=NULL) {
  (length(this$.exclude) > 0)
}, trial=TRUE)


setMethodS3("setExcludedSpots", "MicroarrayData", function(this, excludes, slides=NULL, add=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (any(excludes < 1 || excludes > nbrOfSpots(this)))
    throw("Argument 'excludes' is out of range.");

  slides <- validateArgumentSlides(this, slides=slides);
  
  mat <- matrix(FALSE, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this));
  mat[exludes, slides] <- TRUE;

  if (!hasExcludedSpots(this) || !add) {
    exclude <- mat;
  } else {
    exclude <- unwhich(this$.exclude, dim=dim(mat));
    exclude <- exclude & mat;
  }
  this$.exclude <- which(exlude);
  invisible(this);
}, trial=TRUE)






###########################################################################/**
# @RdocMethod getFieldNames
#
# @title "Gets the names of all fields containing spot specific data"
#
# @synopsis
#
# \description{
#  @get "title". 
#  This is a read-only method. Thus, there is no corresponding set-method.
#  It is used mainly for internal purposes.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFieldNames", "MicroarrayData", function(this) {
  this$.fieldNames;
}, protected=TRUE)







###########################################################################/**
# @RdocMethod getSlideNames
#
# @title "Gets the names of the slides"
#
# @synopsis
#
# \description{
#  @get "title". 
# }
#
# \arguments{
#  \item{slides}{A @vector of slides whose names should be returned. 
#     If @NULL, all slides are considered.}
# }
#
# \value{
#   Returns a @vector of character strings or @NULL, if no slide names are
#   set.
# }
#
# @author
#
# \seealso{
#   @seemethod "setSlideNames".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSlideNames", "MicroarrayData", function(this, slides=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);

  field <- getFieldNames(this)[1];
  
  # Get the current slide names
  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
  colnames <- colnames(this[[field]]);
  colnames[slides];
})


###########################################################################/**
# @RdocMethod setSlideNames
#
# @title "Sets the names of the slides"
#
# @synopsis
#
# \description{
#  @get "title". 
# }
#
# \arguments{
#  \item{names}{A @vector of @character strings for the new slide names.}
#  \item{slides}{A @vector (of the same length as \code{names}) of slides
#     whose names should be set. If @NULL, all slides are considered.}
# }
#
# \value{
#   Returns a @vector of the old slide names that were replaced.
# }
#
# @author
#
# \seealso{
#   @seemethod "getSlideNames".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setSlideNames", "MicroarrayData", function(this, names, slides=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(names)) return();
  if (length(names) != length(slides))
    throw("The number of specified names does not match the number of (specified) slides.");
  
  fields <- getFieldNames(this);
  field <- fields[1];

  # Get the current slide names
  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
  colnames <- colnames(this[[field]]);

  # Remember the old slide names
  oldNames <- colnames[slides];

  # Update the current slide names
  colnames[slides] <- names;
  # No setX() methods nor any static fields should be accessed!
  attr(this, ".memberAccessorOrder") <- c(2,3,5);
  for (field in fields)
    colnames(this[[field]]) <- colnames;

  invisible(oldNames);
})


setMethodS3("getSlideName", "MicroarrayData", function(this, ...) {
  getSlideNames(this, ...);
}, deprecated=TRUE)

setMethodS3("setSlideName", "MicroarrayData", function(this, ...) {
  setSlideNames(this, ...);
}, deprecated=TRUE)




#########################################################################/**
# @RdocMethod getLabel
#
# @title "Gets the label of one field"
#
# @synopsis
#
# \arguments{
#  \item{field}{A field name to get the label for.}
# }
#
# \description{
#   @get "title".
#   Field labels are for instance used as default value is many plot
#   functions. 
#   The labels can be set with the @seemethod "setLabel" method.
# }
#
# \value{
#   Returns the label. If the field was not found, the field is returned.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   print(getLabel(raw, "R"))
# }
#
# @author
#
# \seealso{
#  @seemethod "setLabel".
#  See @see "grDevices::plotmath" how expressions can be used to annotate
#  graphs.
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getLabel", "MicroarrayData", function(this, field) {
  label <- this$.labels[[field]];
  if (is.null(label)) field else label;
})



#########################################################################/**
# @RdocMethod setLabel
#
# @title "Sets the label of one field"
#
# @synopsis
#
# \arguments{
#  \item{field}{A field name to which the label should be assigned.}
#  \item{label}{The new label.}
# }
#
# \description{
#   @get "title" (only). If the field does not exists the
#   label will stored anyway under a virtual field, which could be useful
#   to store different labels.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   setLabel(raw, "R", "red")
#   print(getLabel(raw, "R"))
# }
#
# @author
#
# \seealso{
#   @seemethod "getLabel".
#   See @see "grDevices::plotmath" how expressions can be used to annotate
#   graphs.
#   @seeclass
# }
#*/#########################################################################
setMethodS3("setLabel", "MicroarrayData", function(this, field, value) {
  labels <- this$.labels;
  if (is.null(labels)) labels <- list();  # Should never happens!
  labels[[field]] <- value;
  this$.labels <- labels;
  invisible(this);
})


#########################################################################/**
# @RdocMethod getFlag
#
# @title "Gets the indices of the spots flagged with a specific flag"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{flag}{The name of the flag.}
#  \item{slide}{Slide(s) to be considered. If @NULL all slides are
#   considered.}
# }
#
# \value{
#   Returns a @vector of @integers.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#   setFlag(ma, "|M|>1", abs(ma$M)>1)
#   plot(ma)
#   highlight(ma, "|M|>1", col="red")
# }
#
# @author
#
# \seealso{
#   @seemethod "setFlag",
#   @seemethod "clearFlag",
#   @seemethod "listFlags",
#   @seemethod "getInclude".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getFlag", "MicroarrayData", function(this, flag, slide=NULL) {
  if (missing(flag))
    return(this$Flag);

  flags <- this$.flags;
  if (!is.element(flag, names(flags)))
    throw("No such flag: ", flag);

  dim <- c(nbrOfSpots(this), nbrOfSlides(this));
  res <- matrix(FALSE, nrow=dim[1], ncol=dim[2]);
  idx <- flags[[flag]];
  res[idx] <- TRUE;
  if (is.null(slide)) res else res[,slide];
}, trial=TRUE)


#########################################################################/**
# @RdocMethod setFlag
#
# @title "Sets the flags on some of the spots"
#
# \description{
#   @get "title".
#   See also @seemethod "addFlag" which adds/removes flags to
#   already set flags.
# }
#
# @synopsis
#
# \arguments{
#  \item{flag}{The name of the flag.}
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded.
#   If @NULL no spots are excluded.}
#  \item{slide}{Slide(s) to be considered. If @NULL all slides are
#   considered.}
#  \item{include.op}{If \code{"and"} a spot must be flagged by all the
#   flags to be included, and if the value is \code{"or"} it enough if
#   spot is flagged by only one of the flags.
#   Note that this argument is only effective if the \code{include} 
#   argument(s) are flag names.}
#  \item{exclude.op}{If \code{"and"} a spot must be flagged by all the
#   flags to be excluded, and if the value is \code{"or"} it enough if
#   spot is flagged by only one of the flags.
#   Note that this argument is only effective if the \code{exclude} 
#   argument(s) are flag names.}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#   normalizeWithinSlide(ma, "s")
#   tma <- as.TMAData(ma)
#   setFlag(tma, "Top 5\% M's", top(abs(tma$M), 0.05))
#   setFlag(tma, "Top 5\% T's", top(abs(tma$T), 0.05))
#   setFlag(tma, "Top 5\% M's AND top 5\% T's", c("Top 5\% M's", "Top 5\% T's"))
#   plot(tma)
#   highlight(tma, "Top 5\% M's", col="red")
#   highlight(tma, "Top 5\% T's", col="blue")
#   highlight(tma, "Top 5\% M's AND top 5\% T's", col="purple")
# }
#
# @author
#
# \seealso{
#   @seemethod "addFlag",
#   @seemethod "getFlag",
#   @seemethod "clearFlag",
#   @seemethod "listFlags",
#   @seemethod "getInclude".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("setFlag", "MicroarrayData", function(this, flag, include=NULL, exclude=NULL, slide=NULL, include.op="and", exclude.op="or") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slide <- validateArgumentSlide(this, slide=slide);

  if (is.null(this$.flags)) this$.flags <- list();

  dim <- c(nbrOfSpots(this), nbrOfSlides(this));
  flags <- matrix(FALSE, nrow=dim[1], ncol=dim[2]);

  include <- getInclude(this, include=include, exclude=exclude, include.op=include.op, exclude.op=exclude.op);

  if (!is.null(slide)) {
    if (is.logical(slide)) slide <- witch(slide);
    include[,-slide] <- FALSE;
  }

  this$.flags[[flag]] <- which(flags | include);
  invisible(this);
}, trial=TRUE)


#########################################################################/**
# @RdocMethod addFlag
#
# @title "Flags or unflags the indices of some spots with a specific flag"
#
# \description{
#   @get "title".
#   Note the difference from @seemethod "setFlag" which sets
#   flags regardless of previously flags set, whereas this add/removes to
#   previously set flags.
# }
#
# @synopsis
#
# \arguments{
#  \item{flag}{The name of the flag.}
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded.
#   If @NULL no spots are excluded.}
#  \item{slide}{Slide(s) to be considered. If @NULL all slides are
#   considered.}
#  \item{include.op}{If \code{"and"} a spot must be flagged by all the
#   flags to be included, and if the value is \code{"or"} it enough if
#   spot is flagged by only one of the flags.
#   Note that this argument is only effective if the \code{include} 
#   argument(s) are flag names.}
#  \item{exclude.op}{If \code{"and"} a spot must be flagged by all the
#   flags to be excluded, and if the value is \code{"or"} it enough if
#   spot is flagged by only one of the flags.
#   Note that this argument is only effective if the \code{exclude} 
#   argument(s) are flag names.}
# }
#
# @author
#
# \seealso{
#   @seemethod "getFlag",
#   @seemethod "setFlag",
#   @seemethod "clearFlag",
#   @seemethod "listFlags",
#   @seemethod "getInclude".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("addFlag", "MicroarrayData", function(this, flag, include=NULL, exclude=NULL, slide=NULL, include.op="and", exclude.op="or") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slide <- validateArgumentSlide(this, slide=slide);

  if (is.null(this$.flags)) this$.flags <- list();

  dim <- c(nbrOfSpots(this), nbrOfSlides(this));
  flagged <- this$.flags[[flag]];
  flags <- matrix(FALSE, nrow=dim[1], ncol=dim[2]);
  if (!is.null(flagged))
    flags[flagged] <- TRUE;

  include <- getInclude(this, include=include, exclude=exclude, include.op=include.op, exclude.op=exclude.op);

  if (!is.null(slide)) {
    if (is.logical(slide)) slide <- witch(slide);
    include[,-slide] <- FALSE;
  }

  this$.flags[[flag]] <- which(flags | include);
  invisible(this);
}, trial=TRUE)


#########################################################################/**
# @RdocMethod clearFlag
#
# @title "Clears the spots from the specified flag"
#
# \description{
#   @get "title". The argument \code{flag} can
#   be a regular expression, which means that one can clear more than one
#   flag at the time. All flags are cleared by \code{flag="*"}.
# }
#
# @synopsis
#
# \arguments{
#  \item{flag}{The regular expression of which flags should be cleared.}
# }
#
# @author
#
# \seealso{
#   @seemethod "getFlag",
#   @seemethod "setFlag",
#   @seemethod "listFlags",
#   @seemethod "getInclude".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("clearFlag", "MicroarrayData", function(this, flag) {
  flags <- this$.flags;
  names <- names(flags);
  pos <- regexpr(flag, names);
  idx <- which(pos != -1);
  flags[idx] <- NULL;
  this$.flags <- flags;
  invisible(this);
}, trial=TRUE)


#########################################################################/**
# @RdocMethod listFlags
#
# @title "Lists the names of all known flags"
#
# \description{
#   @get "title" that matches a regular expression.
# }
#
# @synopsis
#
# \arguments{
#  \item{regexpr}{The regular expression of which flags should be cleared.}
# }
#
# @author
#
# \seealso{
#   @seemethod "getFlag",
#   @seemethod "setFlag",
#   @seemethod "clearFlag",
#   @seemethod "getInclude".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("listFlags", "MicroarrayData", function(this, regexpr=NULL) {
  this <- getObject(this);
  flags <- this$.flags;
  names <- names(flags);
  if (!is.null(regexpr)) {
    idx <- regexpr(regexpr, names);
    names <- names[idx != -1];
  }
  names;
}, trial=TRUE)




setMethodS3("getInclude", "MicroarrayData", function(this, include=NULL, exclude=NULL, slides=NULL, include.op="and", exclude.op="or") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);
  
  first <- getFieldNames(this)[1];
  field <- this[[first]];
  field <- pushView(field, getView(this));
  len <- length(field[,1]);
  field <- popView(field);

  dim <- c(len, nbrOfSlides(this));

  # ---------------
  #  I n c l u d e
  # ---------------
  if (include.op == "and")
    inc <- matrix(TRUE, nrow=dim[1], ncol=dim[2])
  else if (include.op == "or")
    inc <- matrix(FALSE, nrow=dim[1], ncol=dim[2])
  else
    throw("Unknown value of argument 'include.op': ", include.op);

  tmp <- NULL;
  if (is.character(include)) {
    for (incl in include) {
      tmp <- this$getFlag(incl);
      if (include.op == "and")
        inc <- inc & tmp 
      else
        inc <- inc | tmp;
    }
    tmp <- NULL;
  } else if (is.matrix(include)) {
    if (nrow(include) != dim[1] || ncol(include) != dim[2])
      throw("When 'include' is a matrix it has to be of the same dimension as the current MicroarrayData object.");
    include <- as.matrix(include[,slides]);
    tmp <- matrix(FALSE, nrow=dim[1], ncol=dim[2]);
    tmp2 <- matrix(FALSE, nrow=dim[1], ncol=length(slides));
    tmp2[include] <- TRUE;
    tmp[,slides] <- tmp2;
    rm(tmp2);
  } else if (is.list(include)) {
    # Since the introduction of Layout$getReplicates() list could now also
    # be used as include/exclude.
    include <- unlist(include);
  } else if (!is.null(include)) {
    tmp <- matrix(FALSE, nrow=dim[1], ncol=dim[2]);
    tmp[include,] <- TRUE
  }
  
  if (!is.null(tmp)) {
    if (include.op == "and")
      inc <- inc & tmp 
    else
      inc <- inc | tmp;
  }


  # ---------------
  #  E x c l u d e
  # ---------------
  if (exclude.op == "and")
    exc <- matrix(TRUE, nrow=dim[1], ncol=dim[2])
  else if (exclude.op == "or")
    exc <- matrix(FALSE, nrow=dim[1], ncol=dim[2])
  else
    throw("Unknown value of argument 'exclude.op': ", exclude.op);

  tmp <- NULL;
  if (is.character(exclude)) {
    for (excl in exclude) {
      tmp <- this$getFlag(excl);
      if (exclude.op == "and")
        exc <- exc & tmp 
      else
        exc <- exc | tmp;
    }
    tmp <- NULL;
  } else if (is.matrix(exclude)) {
    if (nrow(exclude) != dim[1] || ncol(exclude) != dim[2])
      throw("When 'exclude' is a matrix it has to be of the same dimension as the current MicroarrayData object.");
    exclude <- as.matrix(exclude[,slides]);
    tmp <- matrix(FALSE, nrow=dim[1], ncol=dim[2]);
    tmp2 <- matrix(FALSE, nrow=dim[1], ncol=length(slides));
    tmp2[exclude] <- TRUE;
    tmp[,slides] <- tmp2;
    rm(tmp2);
  } else if (is.list(exclude)) {
    # Since the introduction of Layout$getReplicates() list could now also
    # be used as include/exclude.
    exclude <- unlist(exclude);
  } else if (!is.null(exclude)) {
    tmp <- matrix(FALSE, nrow=dim[1], ncol=dim[2]);
    tmp[exclude,] <- TRUE
  }

  if (!is.null(tmp)) {
    if (exclude.op == "and")
      exc <- exc & tmp 
    else
      exc <- exc | tmp;
  }

  as.matrix((inc & !exc)[,slides]);
})


#########################################################################/**
# @RdocMethod as.data.frame
#
# @title "Converts the object to a data frame"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{optional}{If @TRUE, optional fields are also returned.}
# }
#
# \value{
#   Return a @data.frame.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   df <- as.data.frame(raw)
# }
#
# \seealso{
#  @seemethod "extract",
#  @see "base::data.frame".
#   @seeclass
# }
#
# @author
#*/#########################################################################
setMethodS3("as.data.frame", "MicroarrayData", function(x, ..., optional=TRUE) {
  # To please R CMD check...
  this <- x;

  fields <- NULL;
  if (optional == TRUE)
    fields <- c("slide", "spot", "gene", getFieldNames(this));

  extract(this, fields=fields);
})






#########################################################################/**
# @RdocMethod select
#
# @title "Selects rows with certain criteria"
#
# @synopsis
#
# \arguments{
#   \item{...}{The criteria (of fields etc) for returning a row. The indices
#              of all matching rows are returned.}
# }
#
# \description{
#   Given a search criteria, the indices of the rows that fullfills the
#   search criteria will be returned.
# }
#
# \value{Returns a @vector of indices.}
#
# @author
#
# \seealso{
#   @seemethod "extract".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("select", "MicroarrayData", function(this, ...) {
  attach(this);
  on.exit(detach(this));

  if (is.language(...))
    idx <- which(eval(...))
  else
    idx <- which(eval(substitute(...)));

  idx;
})




#########################################################################/**
# @RdocMethod extract
#
# @title "Gets a subset of data specified by fields, slides and/or spot indices"
#
# @synopsis
#
# \description{
#   @get "title".
#   Returns the specified subset of the microarray data as data frame with
#   a first column specifying the slide.
#   If the \code{fields} argument is not specified, all fields are considered.
#   If the \code{slides} argument is not specified, all slides are considered.
#   All data elements returned matches by field (name or index) AND by slide
#   number AND by spot index.
# }
#
# \arguments{
#   \item{fields}{The field names or indices to be returned. If @NULL,
#     all fields are included.}
#   \item{slides}{The slides that should be included. If @NULL,
#     all slides are included.}
#   \item{index}{The spot indices that should be included. If @NULL, 
#     all spots are included.}
#   Any of the above arguments are optional.
# }
#
# \value{Returns a @data.frame.}
#
# \details{
#   This method was formerly named \code{get()}, but since this is a
#   insecure name of a method it was renamed. In addition to the fields
#   as of \code{getFieldNames()}, there are also three special fields:
#   1) \code{slide}, 2) \code{spot} and 3) \code{gene}. These fields are
#   not accessable fields per se, but they are calculated on request.
#   The field \code{gene} is only returned if there are within-slide
#   replicates defined (see class @see "Layout").
#   For instance, \code{as.data.frame()} asks for these fields internally.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   # Get the log ratios (M) for spots 4-7 in the slide 2,3 and 4.
#   extract(ma, c("slide", "M"), slides=2:4, index=4:7)
#   # Gives:
#   #      slide           M
#   #   1      2 -0.41129383
#   #   2      2 -0.44570800
#   #   3      2 -0.16409736
#   #   4      2 -0.22462037
#   #   5      3 -0.51599488
#   #   6      3 -0.30718292
#   #   7      3 -0.54794073
#   #   8      3 -0.42497673
#   #   9      4  0.51019368
#   #   10     4  0.31210389
#   #   11     4  0.08827923
#   #   12     4 -0.12370050
# }
#
# @author
#
# \seealso{
#   @seemethod "as.data.frame".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("extract", "MicroarrayData", function(this, fields=NULL, slides=NULL, index=NULL) {
  virtualFields <- c("slide", "spot", "gene");
  
  if (is.null(fields))
    fields <- getFieldNames(this)
#  else
#    fields <- intersect(fields, c(getFieldNames(this), virtualFields));
          
  # Make sure any fields are retrieved before the get-() and set-() methods.
  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);

  if (is.null(index)) {
    first <- setdiff(fields, virtualFields)[1];
    field <- this[[first]];
    field <- pushView(field, getView(this));
    nbrOfRows <- length(field[,1]);
    field <- popView(field);
    index  <- seq(nbrOfRows);
  } else {
    nbrOfRows  <- length(index);
  }

  slides <- validateArgumentSlides(this, slides=slides);

  nbrOfSlides <- length(slides);
  
  df <- list();

  colnames <- c();
  if ("slide" %in% fields) {
    colnames <- c(colnames, "slide");
    df[["slide"]] <- as.vector(matrix(slides,
                        nrow=nbrOfRows, ncol=nbrOfSlides, byrow=TRUE));
  }

  if ("spot" %in% fields) {
    colnames <- c(colnames, "spot");
    df[["spot"]] <-  as.vector(matrix(index,
                        nrow=nbrOfSlides, ncol=nbrOfRows, byrow=FALSE));
  }

#   if ("gene" %in% fields) {
#     if (hasLayout(this) && hasReplicates(layout <- getLayout(this))) {
# 	colnames <- c(colnames, "gene");
# 	gene <- getGeneIndex(layout, spots=index);
# 	df[["gene"]] <- rep(gene, length.out=length(df[[1]]));
# 	rm(gene);
#     }
#   }

  for (name in setdiff(fields, virtualFields)) {
    colnames <- c(colnames, name);
    field <- this[[name]];
    field <- pushView(field, getView(this));
    df[[name]] <- as.vector(field[index,slides]);
    field <- popView(field);
  }

  df <- as.data.frame(df);
  names(df) <- colnames;
  df;
})





#########################################################################/**
# @RdocMethod nbrOfFields
#
# @title "Gets the number of fields"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfSlides",
#   @seemethod "nbrOfSpots",
#   @seemethod "size".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfFields", "MicroarrayData", function(this) {
  length(getFieldNames(this));
})




#########################################################################/**
# @RdocMethod nbrOfSlides
#
# @title "Gets the number of slides"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfFields",
#   @seemethod "nbrOfSpots",
#   @seemethod "size".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfSlides", "MicroarrayData", function(this) {
  field <- getFieldNames(this)[1];
  if (is.null(field)) return(0);
  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
  value <- this[[field]];
  if (is.null(value))
    0
  else
    ncol(value);
})



#########################################################################/**
# @RdocMethod nbrOfSpots
#
# @title "Gets the number of spots in each of the slides"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfFields",
#   @seemethod "nbrOfSlides",
#   @seemethod "size".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfSpots", "MicroarrayData", function(this) {
  if (hasLayout(this))
    nbrOfSpots(getLayout(this))
  else if (!is.null(field <- getFieldNames(this)[1])) {
    attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
    nrow(this[[field]])
  } else
    0;
})



setMethodS3("nbrOfDataPoints", "MicroarrayData", function(this) {
  if (!is.null(field <- getFieldNames(this)[1])) {
    attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);
    nrow(this[[field]][,])
  } else {
    0;
  }
})


#########################################################################/**
# @RdocMethod size
#
# @title "Gets the number of spots in all slides together"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfFields",
#   @seemethod "nbrOfSlides",
#   @seemethod "nbrOfSpots".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("size", "MicroarrayData", function(this) {
  nbrOfSlides(this) * nbrOfSpots(this);
})




#########################################################################/**
# @RdocMethod hasLayout
#
# @title "Checks if the layout has been specified"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \value{
#   Returns @TRUE if a Layout object has been specified, otherwise @FALSE.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   print(hasLayout(raw))  # TRUE
# }
#
# @author
#*/#########################################################################
setMethodS3("hasLayout", "MicroarrayData", function(this) {
  !is.null(this$layout);
})



#########################################################################/**
# @RdocMethod getLayout
#
# @title "Gets the layout"
#
# @synopsis
#
# \description{
#   Gets the information about the layout structure of the microarray and
#   return it as a \link{Layout} object.
# }
#
# \value{
#   Returns a @see "Layout" object.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   print(getLayout(raw))
# }
#
# @author
#*/#########################################################################
setMethodS3("getLayout", "MicroarrayData", function(this) {
  if (hasLayout(this))
    this$layout
  else
    NULL;
})





#########################################################################/**
# @RdocMethod setLayout
#
# @title "Sets the layout"
#
# \description{
#   Sets the information about the layout structure of the microarray using
#   a \link{Layout} object.
# }
#
# @synopsis
#
# \arguments{
#   \item{layout}{A new \link{Layout} object.}
# }
#
# \examples{
#   SMA$loadData(c("mouse.data", "mouse.setup", "mouse.gnames"))
#   raw <- RawData(mouse.data)
#   setLayout(raw, as.Layout(mouse.setup, id=mouse.gnames))
# }
#
# @author
#*/#########################################################################
setMethodS3("setLayout", "MicroarrayData", function(this, layout) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(layout)) {
    if(!inherits(layout, "Layout"))
      throw("Argument 'layout' must be of class Layout: ", data.class(layout));
  }

  this$layout <- layout;
  invisible(this);
})



setMethodS3("getTreatments", "MicroarrayData", function(this) {
  this$.treatments;
})


setMethodS3("setTreatments", "MicroarrayData", function(this, newTreatments, slides=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);

  len <- length(newTreatments);
  if (length(slides) != length(newTreatments))
    throw("The number of specified 'slides' and the number of treatments do not match.");

  this$.treatments[slides] <- newTreatments;
})


setMethodS3("nbrOfTreatments", "MicroarrayData", function(this) {
  length(unique(this$.treatments))
})




#########################################################################/**
# @RdocMethod getGeneSlideReplicateIndex
#
# @title "Returns the index of the spots in (gene,slide,replicate) indices"
#
# @synopsis
#
# \description{
#   Returns the index of the spots in (spot,slide) indices as
#   (gene,slide,replicate) indices.
# }
#
# \value{
#   Returns a @see "GSRArray" of three dimensions, gene, slide and replicate.
# }
#
# \examples{
#   # Create a raw data object from the preexisting example data in
#   # the sma package.
#   SMA$loadData("mouse.data")
#
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   # Spot 1 & 2, 3 & 4 and so on are replicates of the same gene.
#   setReplicates(layout, "neighboring-pairs")
#
#   raw <- RawData(mouse.data, layout=layout)
#
#   # Get the signal (here by default non-background corrected)
#   ma <- getSignal(raw)
#
#   gsr <- getGeneSlideReplicateIndex(ma)
#
#   #  Replate 1:                Replicate 2:
#   #  ---------------------     ---------------------
#   #  , , 1                     , , 2               
#   #                                                
#   #        [,1] [,2] [,3]            [,1] [,2] [,3]
#   #   [1,]    1   37   73       [1,]    2   38   74
#   #   [2,]    3   39   75       [2,]    4   40   76
#   #   [3,]    5   41   77       [3,]    6   42   78
#   #   [4,]    7   43   79       [4,]    8   44   80
#   #   [5,]    9   45   81       [5,]   10   46   82
#   #   [6,]   11   47   83       [6,]   12   48   84
#   #   [7,]   13   49   85       [7,]   14   50   86
#   #   [8,]   15   51   87       [8,]   16   52   88
#   #   [9,]   17   53   89       [9,]   18   54   90
#   #  [10,]   19   55   91      [10,]   20   56   92
#   #  [11,]   21   57   93      [11,]   22   58   94
#   #  [12,]   23   59   95      [12,]   24   60   96
#   #  [13,]   25   61   97      [13,]   26   62   98
#   #  [14,]   27   63   99      [14,]   28   64  100
#   #  [15,]   29   65  101      [15,]   30   66  102
#   #  [16,]   31   67  103      [16,]   32   68  104
#   #  [17,]   33   69  105      [17,]   34   70  106
#   #  [18,]   35   71  107      [18,]   36   72  108
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/#########################################################################
setMethodS3("getGeneSlideReplicateIndex", "MicroarrayData", function(this) {
  layout <- getLayout(this);
  index <- getGeneReplicateIndex(layout);
  
  gsr <- GSRArray(nbrOfGenes=nrow(index),
                  nbrOfSlides=nbrOfSlides(this),
                  nbrOfReplicates=ncol(index));
  
  for (slide in seq(nbrOfSlides(this)))
    gsr[,slide,] <- index + (slide-1)*nbrOfSpots(this);

  gsr;
}, private=TRUE)




#########################################################################/**
# @RdocMethod setExtra
#
# @title "Sets an extra field"
#
# @synopsis
#
# \arguments{
#   \item{key}{The name under which the extra field should be stored.}
#   \item{value}{The value to be stored.}
# }
#
# \description{
#   Stores optional and extra information of any format.
# }
#
# \examples{
#   raw <- RawData()
#   setExtra(raw, "date", date())
# }
#
# @author
#*/#########################################################################
setMethodS3("setExtra", "MicroarrayData", function(this, key, value) {
  extras <- this$.extras;
  if (is.null(extras)) extras <- list(); # Should never happen!
  extras[[key]] <- value;
  this$.extras <- extras;
  invisible(this);
}, trial=TRUE)


#########################################################################/**
# @RdocMethod getExtra
#
# @title "Gets an extra field"
#
# @synopsis
#
# \arguments{
#   \item{key}{The name under which the extra field is stored.}
# }
#
# \description{
#   Gets an optional and extra field.
# }
#
# \examples{
#   raw <- RawData()
#   setExtra(raw, "date", date())
#   getExtra(raw, "date")
# }
#
# @author
#*/#########################################################################
setMethodS3("getExtra", "MicroarrayData", function(this, key) {
  this$.extras[[key]];
}, trial=TRUE)



setMethodS3("range", "MicroarrayData", function(this, fields, slide=NULL, na.rm=TRUE, inf.rm=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!any(is.element(fields, getFieldNames(this))))
    throw("None of the specified fields are known: ", fields);
  
  X <- extract(this, fields, slide=slide);
  if (inf.rm) X <- X[!is.infinite(X)];
  range(X, na.rm=na.rm);
})



setMethodS3("seq", "MicroarrayData", function(this, ...) {
  seq(length=nbrOfSlides(this));
})




############################################################################
# @RdocMethod getSlidePairs
#
# @title "Get all possible pairs of slides"
#
# \description{
#  @get "title".
# }
#
# \arguments{
#  \item{slides}{A @vector of slides to be included. If @NULL, all slides
#    are considered.}
# }
#
# \value{
#   Returns a @matrix with 2 rows and \eqn{K*(K-1)/2} columns where \eqn{K}
#   is the number of slides to be paired.
#   Column vector \eqn{k} contains the indices of the first and the second
#   slide in the \eqn{k}:th pair.
# }
#
# \examples{
#   # Create a dummy RGData object with K=5 slides
#   X <- matrix(1:10, nrow=2, ncol=5)
#   rg <- RGData(R=X, G=X)
#   pairs <- getSlidePairs(rg)
#   print(pairs)
# }
#
# @author
#
# \seealso{
#  @seeclass
# }
############################################################################
setMethodS3("getSlidePairs", "MicroarrayData", function(this, slides=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);
  nslides <- length(slides);
  if (nslides < 2) {
    throw("No slide pairs available since less than two slides were specified: ", paste(slides, collapse=", "));
  }

  pair <- NULL;
  for (ii in 1:(nslides-1)) {
    for (jj in (ii+1):nslides) {
      pair <- cbind(pair, c(ii,jj));
    }
  }

  apply(pair, MARGIN=2, FUN=function(idx) slides[idx]);
}) # getSlidePairs()




############################################################################
############################################################################
##
##  OBSOLETE?!?
##
############################################################################
############################################################################

setMethodS3("getExtreme", "MicroarrayData", function(this, crit1=0.05, crit2=crit1, what="M", slide=1, include=NULL, exclude=NULL) {

  include <- which(getInclude(this, include, exclude, slide=slide));

  X <- what;
  # To make this generic method a little bit faster, first try to find the
  # data variable 'X' among the fields of the object, i.e. this[[X]].
  # If such a field doesn't exists go an use the get() method instead.
  if (!is.element(X, names(this))) {
    df <- extract(this, fields=X, slide=slide);
    x <- df[,X];
  } else {  
    x <- this[[X]][,slide];
  }
  if (!is.null(include)) x[-include] <- NA;
  
  # Set the upper and lower bounds...
  n <- length(x) - sum(is.na(x));
  if (crit1 >= 1) crit1 <- crit1/n;
  if (crit2 >= 1) crit2 <- crit2/n;
  # Get the extreme spots according to these bounds
  extremes <- (x > quantile.na(x, probs=1-crit2)) | 
              (x < quantile.na(x, probs=crit1));
  which(extremes);
})



setMethodS3("range2", "MicroarrayData", function(this, na.rm=TRUE, inf.rm=TRUE) {
  range <- apply(as.data.frame(this), MARGIN=2, FUN=function(x) {
    if (inf.rm) x <- x[!is.infinite(x)];
    range(x, na.rm=na.rm);
  });
  if (length(range) == 1) range <- c(NA,NA);
  range;
})


setMethodS3("quantile", "MicroarrayData", function(x, probs=seq(0, 1, 0.25), na.rm=TRUE, inf.rm=TRUE, names=TRUE) {
  # To please R CMD check...
  this <- x;

  quantile <- apply(as.data.frame(this), MARGIN=2, FUN=function(x) {
    if (inf.rm) x <- x[!is.infinite(x)];
    quantile(x, probs=probs, na.rm=na.rm, names=names);
  });
  if (length(quantile) == 1) quantile <- rep(NA, length(probs));
  quantile;
})





setMethodS3("getSpotPosition", "MicroarrayData", function(this, slides=NULL, index=NULL) {
  NULL;
})







  
setMethodS3("dataFrameToList", "MicroarrayData", function(this, df, reqFields) {
  df.names <- names(df);

  # Extract or generate the slide index for each data points
  if (is.element("slide", df.names)) {
    slide <- df$slide;
    nbrOfSlides <- diff(range(slide))+1
    dslide <- diff(sort(slide));
    nbrOfSpots <- max(diff(which(dslide != 0)), na.rm=TRUE);
  } else {
    slide <- NULL;
    nbrOfSlides <- NULL;
  }

  # Requires fields...
  fields <- reqFields;

  # Assert that all fields are in the data
  if (!all(is.element(fields, df.names))) {
    throw("Data does not contain all of the required fields: ", fields, collapse=", ");
  }
  
  # Extract or generate the spot index for each data points
  if (is.element("spot", df.names)) {
    spot <- df$spot;
    if (is.null(nbrOfSlides)) {
      # Guess the number of slides.
      dspot <- diff(sort(spot));
      nbrOfSlides <- max(diff(which(dspot != 0)), na.rm=TRUE);
    }
    nbrOfSpots <- max(spot, na.rm=TRUE);
  } else {
    nbrOfSpots <- NULL;
    spot <- NULL;
  }

  # If not found, set the default number of slides.
  if (is.null(nbrOfSlides)) {
    if (inherits(layout, "Layout")) {
      nbrOfSlides <- nbrOfSlides(layout);
    } else {
      warning("Could not find information about the number of slides; will assume there is only one slide.")
      nbrOfSlides <- 1;
    }
  }
  
  # If not found, set the default number of spots.
  if (is.null(nbrOfSpots)) {
    warning("Could not find information about the number of spots; will assume there is as many spots as there is rows in the data divided by the number of slides.")
    nbrOfSpots <- ceiling(nrow(df) / nbrOfSlides);
  }

  # If neither of 'slide' nor 'spot' is specified generated them and
  # construct a vector of indices and from them get the order of these
  # indices.
  if (!is.null(slide) || !is.null(spot)) {
    if (is.null(spot)) {
      spot <- seq(nbrOfSpots);
      spot <- rep(spot, length.out=ncol(df));
    }
    
    if (is.null(slide)) {
      slide <- seq(nbrOfSlides);
      slide <- rep(slide, length.out=ncol(df));
    }
    
    idx <- (slide-1)*nbrOfSpots + spot;
    rm(slide); rm(spot);  # Not of interest anymore.
  } else {
    idx <- NULL;
  }

  # Retrieve the actual field values.
  res <- list();
  for (field in fields) {
    if (!is.element(field, df.names))
      throw("Data does not contain required field: ", field);
  
    if (is.null(idx)) {
      res[[field]] <- matrix(df[[field]], ncol=nbrOfSlides);
    } else {
      res[[field]] <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfSlides);
      res[[field]][idx] <- df[[field]];
    }
  } # for (field ...)
  rm(idx);
  rm(df);

  res;
}, protected=TRUE, static=TRUE, trial=TRUE);





setMethodS3("getBlank", "MicroarrayData", function(this, slides=NULL, include=NULL, blanks="blank|empty", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # 
  layout <- getLayout(this);
  blanks <- getBlank(layout, blanks=blanks)[include];
  matrix(blanks, nrow=length(include), ncol=length(slides));
})




# Examples:
# standardize <- FUN=function(x) {
#   ok <- !is.na(x);
#   c <- median(x[ok]); s <- median(abs(x[ok]-c));
#   (x-c)/s;
# }
#
# v <- applyGroupwise(ma, "M", pd, FUN=standardize, unlist=TRUE);
#
# v <- applyGroupwise(ma, "M", pd, FUN=median, na.rm=TRUE, unlist=TRUE);

setMethodS3("applyGroupwise", "MicroarrayData", function(this, field, layoutGroups, FUN, ..., groups=NULL, slides=NULL, unlist=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);

  spots <- getSpots(layoutGroups, groups=groups);
  data <- this[[field]][,slides];
  l <- list();
  for (k in 1:length(spots)) {
    x <- data[spots[[k]]];
    l[[k]] <- FUN(x, ...);
  }
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})



setMethodS3("applyPrintdipwise", "MicroarrayData", function(this, field, ...) {
  layoutGroups <- getPrintdipGroups(getLayout(this));
  applyGroupwise(this, field=field, layoutGroups=layoutGroups, ...);
})

setMethodS3("applyPrinttipwise", "MicroarrayData", function(this, field, ...) {
  layoutGroups <- getPrinttipGroups(getLayout(this));
  applyGroupwise(this, field=field, layoutGroups=layoutGroups, ...);
})

setMethodS3("applyGenewise", "MicroarrayData", function(this, field, ...) {
  layoutGroups <- getGeneGroups(getLayout(this));
  applyGroupwise(this, field=field, layoutGroups=layoutGroups, ...);
})

setMethodS3("applyPlatewise", "MicroarrayData", function(this, field, ...) {
  layoutGroups <- getPlateGroups(getLayout(this));
  applyGroupwise(this, field=field, layoutGroups=layoutGroups, ...);
})


setMethodS3("getSpotValue", "MicroarrayData", function(this, field, slide=1) {
  # By default, assume that the data is stored spotwise.
  as.matrix(this[[field]])[,slide];
})


setMethodS3("clearCache", "MicroarrayData", function(this, fields=NULL) {
  if (is.null(fields))
    fields <- getFields(this, private=TRUE);
  cachedFields <- fields[grep("^[.]cached[.]", fields)];
  for (field in cachedFields) {
    this[[field]] <- NULL;
  }
  gc();
})

setMethodS3("getCache", "MicroarrayData", function(this, name, force=FALSE) {
  if (force == TRUE)
    return(NULL);
  cacheName <- paste(".cached.", name, sep="");
  value <- this[[cacheName]];
  if (!is.null(value) && !inherits(value, "Cache"))
    warning(paste("Excepted a Cache object (", name, "): ", data.class(value), sep=""));
  value;
})

setMethodS3("setCache", "MicroarrayData", function(this, name, value) {
  cacheName <- paste(".cached.", name, sep="");
  if (!is.null(value) && !inherits(value, "Cache")) {
    class(value) <- unique(c("Cache", class(value)));
  }
  this[[cacheName]] <- value;
})




############################################################################
# HISTORY:
# 2006-05-22
# o seq() now uses seq(length=...) internally, to deal with zero lengths.
# 2006-02-08
# o Rd bug fix: Replaced section 'values' with 'value'.
# o Rd bug fix: Replaced section 'argument' with 'arguments'.
# 2005-03-02
# o Cleaned up the Rdoc comments and harmonized the methods with the
#   generic functions.
# 2005-02-09
# o Updated equal() to return an attribute 'reason' to explain the reason
#   for the detected inequality.
# o Now also append() and equal() support signal and probe weights.
# o Added getChannelNames() to MicroarrayData. Already in GenePixData.
# 2005-02-08
# o keep- and removeSlides() and keep- and removeSpots() now support
#   signal and probe weights.
# 2004-02-17
# o Added keepSpots() and removeSpots() both dated 2003-11-10.
# 2003-12-31
# o Added argument 'fields' to clearCache().
# 2003-12-30
# o Added trial isSaturated() and setSaturated(). Should we use Inf values
#   to indicate saturation instead?!? What's the pro's and con's for that?
#   For instance, Inf values won't be displayed in plots, which means that
#   they have to be replaced by a maximum value all the time. On the other
#   side, Inf's can easily be save and reread by R. If we keep a separate
#   field to keep track of saturated spots we have to save that too. 
#   However, with Inf's we somehow have to specify for each channel what
#   the maximum value is (or assume it is equal to max(65535, max value))
#   when saving and rereading from file.
# o getSlideNames() replaced getSlideName(), which is made deprecated. Same
#   for the set-method. Added Rdoc comments too.
# o Added removeSlides() and keepSlides() with Rdoc comments.
# o Added Rdoc comments for getFieldNames() and summary().
# 2003-10-13
# o BUG FIX: setCache() gave an error is cache value was NULL.
# 2003-09-26
# o Added getCache(), setCache(), clearCache().
# 2003-07-29
# o Updated append() to set the Layout object if such an object does *not* 
#   already exist, but the appended MicroarrayData object has one. 
# o Added an example to append() showing how one easily can read multiple
#   GPR files one by one and convert each to MAData objects that are 
#   appended together. This makes it possible to read 10-20 times more
#   data files.
# 2003-05-04
# o BUG FIX: The example of getGeneReplicateIndex() contained syntax errors.
# 2002-12-24
# o BUG FIX: setSlideName() failed if there existed a virtual field with
#   the same name.
# 2002-12-01
# o range() was duplicated. Renamed the obsolete version to range2().
# 2002-11-28
# o Added summary().
# 2002-11-27
# o BUG FIX: If both MicroarrayData objects did not have the Layout object
#   set, append() would throw an Exception. Now it is possible to append
#   to object with no Layout set.
# 2002-11-12
# o Updated getField() by changing argument 'view' to 'viewFunction', which
#   can be NULL, as.vector(), as.SpotSlideArray(), as.GeneSlideArray() etc.
#   BTW: Should getField() be made obsolete. Is it actually needed?
# 2002-11-06
# o Added nbrOfDataPoints()
# 2002-10-31
# o Added getSpotValue().
# 2002-10-09
# o BUG FIX: as.character() would create a vector of string with the length
#   equal to the number of fields, one string for each field.
# 2002-09-30
# o Made (all) MicroarrayData classes Morphable.
# 2002-08-06
# o If equals() returns FALSE, FALSE got an attribute "firstUnequalField"
#   showing which field that was unequal. Note that more fields might
#   differ.
# 2002-07-14
# * Updated getInclude() to support 'include' and 'exclude' as matrices too.
#   It also is guaranteed to return a matrix (even if has one column).
# 2002-07-02
# * Updated extract() to return the true field names. Before there where
#   some automatic 'safe data frame names' translation.
# 2002-06-24
# * Made set() protected.
# 2002-05-06
# * All normalization methods are now in MicroarrayData.NORM.R.
# * Added a assertion that argument layout in the constructor is of class
#   Layout.
# 2002-05-04
# * Added applyGroupwise() and the applyXXXwise()'s.
# 2002-05-03
# * Added getBlanks().
# 2002-04-28
# * BUG FIX: adjustBiasScale...() did not scale correctly if bias was NOT
#   adjusted!
# 2002-04-21
# * Extracted log functions to MicroarrayData.LOG.R.
# * Extracted I/O functions to MicroarrayData.IO.R.
# * Extracted plot functions to MicroarrayData.PLOT.R.
# * Made getColors() generate grayscale colors by default. Before the method
#   was declared abstract.
# * Added trial version of normalizeGenewise().
# * Replaced some of the throw()'s with throw()'s.
# * Added getGeneReplicateIndex() and getGeneSlideReplicateIndex()
#   to MicroarrayData for fast access to the (gene, slide, replicate)
#   indices.
# 2002-04-20
# * Added trial versions of the static functions dataFrameToList() and
#   readToList(). These can support the read()  functions in the subclasses.
# * hist() now also excludes Inf's in addition to NA's. It turned out that
#   Lei Jiang's data, who reported the bug,  contained *one* M value that
#   was Inf. The result was that pretty(), which is called by hist(), gave
#   an error like "NA/NaN/Inf in foreign function call (arg 2)".
# * Added reference to 'plot' also in addition to 'par' in plot-functions.
#   This was done on a question how to set the limits on the axis.
# 2002-04-12
# * Updated the Rdoc example for plot() so it gives example on how to plot
#   a certain slide.
# 2002-04-06
# * Added support for multiple fields in normalizePrintorder() and
#   normalizeSpatial().
# * Renamed argument 'what' to 'field' in plotPrintorder().
# * Added argument 'breakpoints' to plotPrintorder().
# * Added normalizePrintorder().
# 2002-04-05
# * Added normalizeSpatial(). It is "smart", because it tries to get the
#   physical positions of the spots, but such information is not available
#   the positions according to the Layout object is used.
# 2002-04-03
# * write() now supports File objects.
# 2002-03-29
# * Updated the Rdoc's so any references to old get() are now to extract().
# 2002-03-25
# * Added static method createColors().
# 2002-03-24
# * Updated the Rdoc example for text().
# * BUG FIX: text() did not work correctly if argument 'labels' where not
#   explicitly set. The reason for this was that getInclude() was changed
#   to return a vector of TRUE and FALSE. Had to use which().
# * Changed this$getInclude(...) to getInclude(this, ...).
# 2002-03-13
# * BUG FIX: Forgot the sep="\t" in write().
# * Added some more Rdoc comments to write().
# 2002-03-10
# * BUG FIX: Argument 'slides' in write() was mistakenly named 'slide'.
# * write() in MicorarrayData will now by default save the object as the
#   data frame returned by as.data.frame(), i.e. if write() is not
#   overridden by a subclass e.g. GenePixData, which then saves in a
#   different format.
# * BUG FIX: extract() in MicroarrayData would neglect the special fields
#   "slide", "spot" and "gene". This automatically also fixed the fact that
#   as.data.frame() would not create these fields.
# 2002-02-26
# * Added the virtual fields "slide", "spot", "gene" to extract.
# 2002-02-25
# * Added read(), readAll(), write(), append().
# * Modified this class to support GenePixData etc directly without making
#   use of a ResultsData class.
# 2002-01-24
# * Renamed method get() to extract(). get() is not safe!
# * Rewritten to make use of setClassS3 and setMethodS3.
# 2002-01-19
# * Added getSlideName() and setSlideName(). Used in automatic labelling of
#   plots. See the putSlide() method.
# 2002-01-17
# * Added argument 'new=TRUE' to printReplicates. With new==FALSE it is
#   possible to plot another sequence of replicates in the same plot. This
#   is useful for instance when you want to look at the effect before and
#   after normalization.
# * BUG FIX: putGene() crashed if 'id' or 'name' was "auto"; instead of 
#   doing "if(id && name) ..." it is safer/better to do 
#   "if (id == TRUE && name == TRUE) ...".
# * BUG FIX: plotReplicates() didn't work if there where no within-slide
#   replicates. This is now corrected.
# 2002-01-15
# * Added putSlide().
# * Added seq() to simplify multiplots.
# * Added argument 'adjustMargins=TRUE' to subplots().
# * Added putTimestamp, putDataset, and putAuthor to all plot 
#   methods.
# 2001-11-18
# * Added getField() and getFields().
# 2001-08-11
# * Added plotPrintorder.
# 2001-08-08
# * Added the first core functionality of has-, get- setExcludedSpots.
#   By using which and unwhich I hope it could also be memory efficient.
#   I do not want to just set the values to NA to exclude, because then
#   one can not unexclude. Instead I want to flag some spots to be excluded.
# 2001-08-06
# * BUG FIX: plotSpatial set the plot history to "plotXY" instead of
#   "plotSpatial". This bug was probably a cut-and-paste misstake.
# * BUG FIX: plotXY and plotSpatial generated the wrong colors for all
#   cases where slide > 1. It turned out do be bug in getColors.MAData.
# 2001-08-03
# * Added getFieldNames(), setFieldNames(), renameField().
# 2001-07-31
# * Added hasLayout() and made setLayout() assert the argument layout.
# * Updated nbrOfSpots() to either ask the Layout of use 
#   get(fields=1, slides=1) to figure out the number of spots. This
#   method is a little bit inefficient so it should be overloaded by
#   subclasses.
# 2001-07-24
# * lowessCurve() now returns the lowess line.
# 2001-07-18
# * Bug fix: Forgot the 'labels' argument in call to text() in plotXY.
# 2001-07-17
# * TODO: Have to decide if spot indices should also run across slides, i.e.
#   the unique indices should continue counting on the next slide etc. This
#   is a complicated issue and somehow the though has to be digested before
#   making a decision. Now, I think some function are a little bit ambigous
#   in the use of slide, include and exclude arguments. It is pretty obvious
#   though that the include and exclude should be done before applying the
#   slide argument.
# * Added some Rdoc comments.
# * Updated the plot() method internally; now less if and then statements.
# * Removed the .lastPlot field. Making use of Device$setPlotParameters
#   instead.
# 2001-07-15
# * getIncludes() now also accepts lists in include/exclude arguments.
# * From now on the cex, col, pch etc arguments to the plot functions are
#   not subject to include and exclude. I made this decision since most
#   often you want for instance highlight four spots with four different
#   colors and nothing else. This was much harder to do before.
# 2001-07-12
# * Updated the pin.lty in plotXY to be the same as the one in Layout$put().
# * Bug fix: Trying to load a gpr data set with layout 8x4x15x16, the
#   pin.lty <- rep(...) function gave an error. used ngrid.c instead of
#   ngrid.r. The data I've tried this far have had ngrid.c == ngrid.r!
# 2001-07-11
# * Made .layout public, i.e. renamed it to layout.
# * Updated some of the Rdoc comments.
# 2001-07-09
# * Totally removed the use of image() in the plotSpatial() method. image()
#   had a "uncontrolable" color method and thanks to a highly improved
#   getPosition.Layout speed.
# 2001-07-06
# * Added addFlag(). Updated clearFlag() to work with regexpr's too.
# * Renamed plotYvsX to plotXY.
# 2001-07-05
# * Made include in the plot functions be operating on flags which have been
#   set by flag(). Removed all exlude from the plot functions.
# * Made plotYvsX() and plotSpatial() more generic and moved both of them to
#   this class. Also moved highlight(), text(), plot() to this class.
# 2001-07-04
# * Added the .extra field w/ the setExtra(), getExtra() methods.
# 2001-07-01
# * Removed plotSpatial(); now MicroarrayData is totally plot free.
# * Added getLabel() and setLabel(). Really useful!
# * Removed the rename() method since the internal field names should never
#   be modified.
# * Generic get() and set() seems to works fine. Added abstract setField().
# 2001-06-29
# * Created from old ResultsData.
############################################################################
