###########################################################################/**
# @RdocClass AffymetrixDataFile
#
# @title "The abstract AffymetrixDataFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixDataFile object represents a single Affymetrix data file,
#  e.g. an Affymetrix CEL file or an Affymetrix Probe Data (APD) file.
#  Note that this class is abstract and can not be instanciated, but
#  instead you have to use one of the subclasses or the 
#  @seemethod "fromFile" method.
# }
# 
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the data file.}
#   \item{path}{An optional path to the data data file.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
#*/###########################################################################
setConstructorS3("AffymetrixDataFile", function(...) {
  extend(AffymetrixFile(...), "AffymetrixDataFile",
    .apdMap = NULL,
    .doTransform = FALSE,
    .transforms = NULL
  )
}, abstract=TRUE)



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the Affymetrix data file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "AffymetrixDataFile", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Name: ", getName(this), ".", sep="");
  s <- paste(s, " Chip type: ", getChipType(this), ".", sep="");
  if (hasTransforms(this)) {
    s <- paste(s, " Has tranforms: TRUE.", sep="");
    s <- paste(s, " Tranform signals: ", isTransforming(this), ".", sep="");
  }
  if (!is.null(getApdMap(this))) {
    s <- paste(s, " Using reading map: ", getApdMap(this), sep="");
  }
  s <- sprintf("%s Pathname: %s (%.2fMb).", s, getPathname(this), 
                                                   getFileSize(this)/1024^2);
  s <- sprintf("%s RAM: %.2fMb.", s, objectSize(this)/1024^2);
  s;
})




###########################################################################/**
# @RdocMethod readIntensities
#
# @title "Reads probe intensities from a data file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("readIntensities", "AffymetrixDataFile", abstract=TRUE, protected=TRUE);




###########################################################################/**
# @RdocMethod getProbeIntensities
#
# @title "Gets the probe intensities for a set of probes"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{probes}{Am index @vector specifying which probe to be queried. 
#    If @NULL, all probes are considered.}
#  \item{...}{Not used.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{.readMap}{(Internal) A read map}
#  \item{.doTransform}{(Internal) Transform signals or not?}
# }
#
# \value{
#   Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbeIntensities", "AffymetrixDataFile", function(this, probes=NULL, ..., verbose=FALSE, .readMap=getReadMap(this), .doTransform=isTransforming(this)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probes':
  nbrOfProbes <- nbrOfProbes(this);
  if (is.null(probes)) { 
  } else {
    probes <- Arguments$getIndices(probes, range=c(1,nbrOfProbes));
    nbrOfProbes <- length(probes);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Reading probe signals");

  y <- readIntensities(this, probes=probes, readMap=.readMap);   # Abstract
  verbose && exit(verbose);

  if (.doTransform)
    y <- transformProbeSignals(this, y, verbose=verbose);

  y;
})


###########################################################################/**
# @RdocMethod "["
#
# @title "Gets the probe intensities for a set of probes"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{i}{An index @vector specifying which probe to be queried. 
#    If @NULL, all probes are considered.}
#  \item{drop}{If @TRUE and the length of the return vector is less or
#    equal to one, then its dimensions are dropped, otherwise not.}
# }
#
# \value{
#   Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getProbeIntensities".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("[", "AffymetrixDataFile", function(this, i=NULL, drop=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'i':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  y <- getProbeIntensities(this, probes=i);
  if (length(y) < 1 && drop)
    dim(y) <- NULL;
  y;
})




###########################################################################/**
# @RdocMethod getProbesetIntensities
#
# @title "Gets probe signals for a subset of probesets in one array"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{probesets}{An @integer index @vector specifying probesets to be read. 
#    If @NULL, all probesets are read.}
#  \item{doTransform}{If @TRUE, probe signals are transformed by the
#    internal transformation functions, otherwise not.}
#  \item{...}{Arguments passed to the low-level function for read probesets, 
#    e.g. @see "affxparser::readCelUnits" or @see "aroma.apd::readApdUnits".}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a named @list structure.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProbesetIntensities", "AffymetrixDataFile", function(this, probesets, doTransform=isTransforming(this), ..., stratifyBy=c("nothing", "pmmm", "pm", "mm"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'doTransform':
  doTransform <- Arguments$getLogical(doTransform);

  # Argument 'stratifyBy':
  stratifyBy <- match.arg(stratifyBy);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Transform signals?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (doTransform && hasTransforms(this)) {
     # Building transforms functions
     verbose && enter(verbose, "Building transforms functions");
     transforms <- function(y, ...) {
       transformProbeSignals(this, y=y);
     }
     verbose && exit(verbose);
  } else {
     transforms <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  readUnits(this, units=probesets, ..., stratifyBy=stratifyBy, transforms=transforms, verbose=verbose);
})



############################################################################
# HISTORY:
# 2006-08-11
# o Now AffymetrixDataFile inherits from AffymetrixFile.
# 2006-05-22
# o Added more Rdoc comments.
# 2006-05-15
# o Renamed getMap() to getReadMap() and same for the other map functions.
#   Removed mapOn(), mapOff() etc.
# 2006-04-09
# o Added abstract method getMapType().
# 2006-04-05
# o BUG FIX:  fitQuantileNormFcn() returned a function that when called with
#   x values forcing extrapolation, error "Error in Recall(object, xrange) :
#   couldn't find function "predict.smooth.spline.fit" would be thrown.
#   This is because you cannot do predict(sp$fit, ...) but only 
#   predict(sp, ...).  Why I don't really know; probably about namespaces.
# 2006-03-23
# o Added findCdf().
# 2006-03-18
# o Added argument 'subset' to fitQuantileNormFcn().
# 2006-03-04
# o Remapping of indices works only with integers, which means that we can
#   only index 256^4 = 4.2*10^9 cells.
# o Added get- and setMap() and hasMap().
# o Now getProbeSignals() and readIntensities() accepts argument 'map'.
# 2006-03-03
# o When creating transformation function in, say, fitQuantileNormFcn(), it
#   is important to create an empty environment for the function otherwise
#   all arguments in the calling function is included too.
# o Added writeApd().  For now it can only write 'intensities'.
# 2006-03-02
# o Created.  This will allow us to define Affymetrix data sets of any kind,
#   for instance based on CEL files or APD files (or even mixed).
############################################################################
