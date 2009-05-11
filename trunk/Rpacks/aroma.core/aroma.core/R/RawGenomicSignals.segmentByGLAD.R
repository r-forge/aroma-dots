###########################################################################/**
# @set "class=RawGenomicSignals"
# @RdocMethod segmentByGLAD
#
# @title "Segment copy numbers using the GLAD method"
#
# \description{
#  @get "title" of the \pkg{GLAD} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to the segmentation function.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the fit object.
# }
# 
# \details{
#   Internally @seemethod "GLAD::glad" is used to segment the signals.
#   This segmentation method support weighted segmentation.
# }
#
# @author
#
# \seealso{
#   Internally the segmentation function 
#   @seemethod "GLAD::glad" is used.
#   @seeclass
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("segmentByGLAD", "RawGenomicSignals", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Segmenting");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving the fit function");
  pkgName <- "GLAD";
  # Assert that package is installed
  isPackageInstalled(pkgName) || throw("Package is not installed: ", pkgName);
  pkg <- packageDescription(pkgName);
  pkgVer <- pkg$Version;
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer);

  methodName <- "glad";
  verbose && cat(verbose, "Method: ", methodName);
  verbose && cat(verbose, "Package: ", pkgDetails);

  # We need to load package
  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);

  # Get the fit function for the segmentation method
  envir <- as.environment(sprintf("package:%s", pkgName));
  fitFcn <- get(methodName, mode="function", envir=envir);
  verbose && str(verbose, "Function: ", fitFcn);
  formals <- formals(fitFcn);
  verbose && cat(verbose, "Formals:");
  verbose && str(verbose, formals);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data of interest");
  data <- extractDataForSegmentation(this, ..., verbose=less(verbose, 5));
  verbose && str(verbose, data);
  verbose && exit(verbose);

  sampleName <- attr(data, "sampleName");
  chromosome <- data$chromosome[1];
  nbrOfLoci <- nrow(data);
  hasWeights <- !is.null(data$w);

  verbose && cat(verbose, "Sample name: ", sampleName);
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasWeights) {
    # Verify that weights are supported (not yet)
    if (!is.element("weights", names(formals))) {
      hasWeights <- FALSE;
      msg <- paste("Weights detected but ignored, because the available segmentation function ('", methodName, "()') does not support weights. Check with a more recent version of the package: ", pkgDetails, sep="");
      verbose && cat(verbose, "WARNING: ", msg);
      warning(msg);
    }
  } # if (hasWeights)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure");
  cnData <- data.frame(
    LogRatio=data$y, 
    PosOrder=1:nbrOfLoci, 
    Chromosome=data$chromosome,
    PosBase=data$x
    # Add (chipType, units) identifiers to be able to backtrack SNP IDs etc.
#    chipType=as.factor(chipType),
#    units=units,
#    sdTheta=data$sdTheta
  );
  verbose && str(verbose, cnData);
  cnData <- GLAD::as.profileCGH(cnData);
  verbose && str(verbose, cnData);
  verbose && exit(verbose);

  args <- list(cnData);

  if (hasWeights) {
    fitArgs <- list(weights=data$w);
    verbose && cat(verbose, "Additional segmentation arguments:");
    keep <- (names(fitArgs) %in% names(formals));
    fitArgs <- fitArgs[keep];
    verbose && str(verbose, fitArgs);
    args <- c(args, fitArgs);
    rm(fitArgs);
  }

  userArgs <- list(...);
  if (length(userArgs) > 0) {
    verbose && cat(verbose, "User and segmentation arguments:");
    verbose && str(verbose, userArgs);
    verbose && cat(verbose, "Kept user arguments:");
    keep <- (names(userArgs) %in% names(formals));
    userArgs <- userArgs[keep];
    verbose && str(verbose, userArgs);
    args <- c(args, userArgs);
    rm(userArgs);
  }

  verbose && cat(verbose, "Final arguments:");
  args <- c(args, list(verbose=as.logical(verbose)));
  verbose && str(verbose, args);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calling segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Calling %s() of %s", methodName, pkgName));
  # In case the method writes to stdout, we capture it
  stdout <- capture.output({
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    t <- system.time({
      fit <- do.call(methodName, args);
    });
    attr(fit, "processingTime") <- t;
  });

  verbose && cat(verbose, "Captured output that was sent to stdout:");
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && cat(verbose, "Fitting time (in seconds):");
  verbose && print(verbose, t);

  verbose && cat(verbose, "Fitting time per 1000 loci (in seconds):");
  verbose && print(verbose, 1000*t/nbrOfLoci);

  verbose && cat(verbose, "Results object:");
  verbose && str(verbose, fit);

  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;  
}, private=TRUE) # segmentByCBS()


############################################################################
# HISTORY:
# 2009-05-10
# o Created.
############################################################################
