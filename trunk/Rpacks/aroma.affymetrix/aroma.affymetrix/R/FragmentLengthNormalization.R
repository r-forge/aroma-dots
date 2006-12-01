###########################################################################/**
# @RdocClass FragmentLengthNormalization
#
# @title "The FragmentLengthNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for PCR 
#  fragment length effects on total copy-number chip-effect estimates.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{ces}{A @see "CnChipEffectSet".}
#   \item{targetFunction}{A @function.  The target function to which all arrays
#     should be normalized to.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.
#     If @NULL, all are considered.}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \examples{\dontrun{
# }}
#
# @author
#*/###########################################################################
setConstructorS3("FragmentLengthNormalization", function(ces=NULL, targetFunction=NULL, subsetToFit=NULL, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ces':
  if (!is.null(ces)) {
    if (!inherits(ces, "CnChipEffectSet"))
      throw("Argument 'ces' is not an CnChipEffectSet object: ", class(ces));

    if (ces$combineAlleles != TRUE) {
      throw("Currently only total copy-number chip effects can be normalized, i.e. 'combineAlleles' must be TRUE");
    }

    if (ces$mergeStrands != TRUE) {
      throw("Currently only total copy-number chip effects can be normalized, i.e. 'mergeStrands' must be TRUE");
    }
  }

  if (!is.null(targetFunction)) {
    if (!is.function(targetFunction)) {
      throw("Argument 'targetFunction' is not a function: ", class(targetFunction)[1]);
    }
  }

  this <- extend(Object(), "FragmentLengthNormalization", 
    .tags = tags,
    inputSet = ces,
    "cached:outputSet" = NULL,
    .params = list(
       subsetToFit = subsetToFit,
       .targetFunction = targetFunction
    )
  );

  setTags(this, tags);

  this;
})



setMethodS3("as.character", "FragmentLengthNormalization", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  ces <- getInputSet(this);
  tags <- paste(getTags(ces), collapse=", ");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Number of arrays: %d (%.2fMb)", 
                           nbrOfArrays(ces), getFileSize(ces)/1024^2));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ces))));
  params <- paste(getParametersAsString(this), collapse=", ");
  s <- c(s, sprintf("Algorithm parameters: (%s)", params));
  s <- c(s, sprintf("Output path: %s", getPath(this)));
  s <- c(s, sprintf("Is done: %s", isDone(this)));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output data set"
#
# \description{
#  @get "title", which is the same as the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
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
#*/###########################################################################
setMethodS3("getName", "FragmentLengthNormalization", function(this, ...) {
  ces <- getInputSet(this);
  getName(ces);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this normalizer.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "FragmentLengthNormalization", function(this, ...) {
  tags <- this$.tags;

  ds <- getInputSet(this);
  tags <- c(getTags(ds), tags);

  # Update default tags
  tags[tags == "*"] <- "FLN";

  tags;
})


setMethodS3("setTags", "FragmentLengthNormalization", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }
  
  this$.tags <- tags;
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output data set"
#
# \description{
#  @get "title", which is the name with comma separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
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
#*/###########################################################################
setMethodS3("getFullName", "FragmentLengthNormalization", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "FragmentLengthNormalization", function(this, ...) {
  "normFragmentLength";
})


setMethodS3("getParametersAsString", "FragmentLengthNormalization", function(this, ...) {
  params <- getParameters(this);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, protected=TRUE)

setMethodS3("getParameters", "FragmentLengthNormalization", function(this, ...) {
  this$.params;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the output data set"
#
# \description{
#  @get "title".
#  If non-existing, then the directory is created.
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
#*/###########################################################################
setMethodS3("getPath", "FragmentLengthNormalization", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getInputSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);
  chipType <- gsub("-monocell$", "", chipType);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})


###########################################################################/**
# @RdocMethod getInputSet
#
# @title "Gets the source chip effects"
#
# \description{
#  @get "title" that is to be normalized (or has been normalized).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getInputSet", "FragmentLengthNormalization", function(this, ...) {
  this$inputSet;
})



###########################################################################/**
# @RdocMethod getOutputSet
#
# @title "Gets the normalized data set"
#
# \description{
#  @get "title", if normalized.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet" or @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getOutputSet", "FragmentLengthNormalization", function(this, ..., force=FALSE) {
  outputSet <- this$outputSet;
  if (force || is.null(outputSet)) {
    if (isDone(this)) {
      inputSet <- getInputSet(this);
      clazz <- Class$forName(class(inputSet)[1]);
      outputSet <- clazz$fromFiles(path=getPath(this));

      # Ad hoc for now /HB 2006-12-01
      outputSet$mergeStrands <- inputSet$mergeStrands;
      outputSet$combineAlleles <- inputSet$combineAlleles;

      this$outputSet <- outputSet;
    }
  }
  outputSet;
})


setMethodS3("getOutputFiles", "FragmentLengthNormalization", function(this, ...) {
  outPath <- getPath(this);
  findFiles(pattern="[.](c|C)(e|E)(l|L)$", paths=outPath, firstOnly=FALSE);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod isDone
#
# @title "Checks if the data set is normalized or not"
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
#  Returns @TRUE if the data set is normalized, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isDone", "FragmentLengthNormalization", function(this, ...) {
  pathnames <- getOutputFiles(this);
  if (length(pathnames) == 0)
    return(FALSE);

  ds <- getInputSet(this);  
  if (length(pathnames) != nbrOfArrays(ds)) {
    throw("Number of output CEL files does not match the number of CEL files in the input dataset: ", length(pathnames), " != ", nbrOfArrays(ds));
  }
  
  return(TRUE);
})


setMethodS3("getCdf", "FragmentLengthNormalization", function(this, ...) {
  inputSet <- getInputSet(this);
  getCdf(inputSet);
})


setMethodS3("getSubsetToFit", "FragmentLengthNormalization", function(this, force=FALSE, ...) {
  # Cached?
  units <- this$.units;
  if (!is.null(units) && !force)
    return(units);

  # Identify all SNP units
  cdf <- getCdf(this);
  units <- indexOf(cdf, "^SNP");

  # Keep only those for which we have PCR fragmenth-length information
  si <- getSnpInformation(cdf);
  fl <- getFragmentLengths(si, units=units);
  keep <- is.finite(fl);
  units <- units[keep];

  # Fit to a subset of the units?
  subsetToFit <- this$.params$subsetToFit;
  if (!is.null(subsetToFit)) {
    # A fraction subset?
    if (length(subsetToFit) == 1 && 0 < subsetToFit && subsetToFit < 1) {
      keep <- seq(from=1, to=length(units), length=subsetToFit*length(units));
    } else {
      keep <- which(units %in% subsetToFit);
    }

    # Make sure to keep data points at the tails too
    keep <- c(keep, which.min(fl), which.max(fl));
    keep <- unique(keep);

    # Now filter
    units <- units[keep];
  }

  # Sort units
  units <- sort(units);

  # Assert correctness
  units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));

  # Cache
  this$.units <- units;

  units;
}, protected=TRUE)



setMethodS3("getTargetFunction", "FragmentLengthNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  fcn <- this$.targetFunction;
  if (is.null(fcn) || force) {
    verbose && enter(verbose, "Estimating target prediction function");

    # Get target set
    ces <- getInputSet(this);
    ceR <- getAverageFile(ces, verbose=less(verbose));

    # Get units to fit
    units <- getSubsetToFit(this);

    # Get PCR fragment lengths for these
    cdf <- getCdf(this);
    si <- getSnpInformation(cdf);
    fl <- getFragmentLengths(si, units=units);

    # Get target log2 signals for SNPs
str(111);
str(units);
    yR <- getDataFlat(ceR, units=units, fields="theta", verbose=less(verbose));
str(yR);
    yR <- yR[,"theta"];
    yR <- log2(yR);

    # Fit lowess function
    verbose && enter(verbose, "Fitting target prediction function");
    ok <- (is.finite(fl) & is.finite(yR));
    fit <- lowess(fl[ok], yR[ok]);
    class(fit) <- "lowess";

    # Create target prediction function
    fcn <- function(x, ...) {
      predict(fit, x, ...);  # Dispatched predict.lowess().
    }
    verbose && exit(verbose);

    verbose && exit(verbose);

    this$.targetFunction <- fcn;
  }

  fcn;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized, 
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
setMethodS3("process", "FragmentLengthNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing set for PCR fragment-length effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputSet <- getOutputSet(this);
    return(invisible(outputSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve/calculate the target function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  targetFcn <- getTargetFunction(this, verbose=less(verbose));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input dataset
  ces <- getInputSet(this);

  # Get SNP units
  cdf <- getCdf(ces);
  subsetToUpdate <- indexOf(cdf, "^SNP");

  verbose && enter(verbose, "Retrieving PCR fragment lengths");
  cdf <- getCdf(this);
  si <- getSnpInformation(cdf);

  # Get PCR fragment lengths for the subset to be fitted
  fl <- getFragmentLengths(si, units=subsetToUpdate);
  verbose && exit(verbose);

  # Get subset to fit
  subsetToFit <- getSubsetToFit(this, verbose=less(verbose));

  # Get (and create) the output path
  path <- getPath(this);
  mkdirs(path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  map <- NULL;
  res <- vector("list", nbrOfArrays(ces));
  for (kk in seq(length=nbrOfArrays(ces))) {
    ce <- getFile(ces, kk);
    verbose && enter(verbose, sprintf("Array #%d (%s)", kk, getName(ce)));

    if (is.null(map)) {
      verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
      map <- getCellMap(ce, units=subsetToUpdate, verbose=less(verbose));
      verbose && str(verbose, map);
      verbose && exit(verbose);
    }

    filename <- getFilename(ce);
    pathname <- filePath(path, filename);
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized. Skipping.");
      ceN <- fromFile(ce, pathname);
      # CDF inheritance
      setCdf(ceN, cdf);
      res[[kk]] <- ceN;
      verbose && exit(verbose);
    }

    # Get target log2 signals for all SNPs to be updated
    verbose && enter(verbose, "Getting signals");
    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
    verbose && exit(verbose);

    # Extract the values to fit the normalization function
    verbose && enter(verbose, "Normalizing log2 signals");
    y <- data[,"theta"];
    subset <- match(subsetToFit, subsetToUpdate);
    yN <- normalizeFragmentLength(log2(y), fragmentLengths=fl, 
                             targetFcn=targetFcn, subsetToFit=subset, ...);
    yN <- 2^yN;
    verbose && exit(verbose);

    # Copy CEL file and update the copy
    verbose && enter(verbose, "Copying source CEL file");
    copyCel(from=getPathname(ce), to=pathname);
    verbose && exit(verbose);

    # Defining normalized object
    ceN <- fromFile(ce, pathname);
    # CDF inheritance
    setCdf(ceN, cdf);

    verbose && enter(verbose, "Storing normalized signals");
    data[,"theta"] <- yN;
    updateDataFlat(ceN, data=data, verbose=less(verbose));
    verbose && exit(verbose);

    res[[kk]] <- ceN;

    verbose && exit(verbose);
  }

  # Create the output set (ad hoc for now so that we keep parameter too)
  outputSet <- clone(ces);
  outputSet$files <- res;

  # Update the output dataset
  this$outputSet <- outputSet;

  verbose && exit(verbose);
  
  outputSet;
})

############################################################################
# HISTORY:
# 2006-11-28
# o Created from QuantileNormalizer.R.
############################################################################
