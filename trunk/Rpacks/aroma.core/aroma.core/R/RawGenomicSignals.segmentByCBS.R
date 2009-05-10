###########################################################################/**
# @set "class=RawGenomicSignals"
# @RdocMethod segmentByCBS
#
# @title "Segment copy numbers using the CBS method"
#
# \description{
#  @get "title" of the \pkg{DNAcopy} package.
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
#   Internally @seemethod "DNAcopy::segment" is used to segment the signals.
#   This segmentation method support weighted segmentation.
# }
#
# @examples "../incl/RawGenomicSignals.SEG.Rex"
#
# @author
#
# \seealso{
#   Internally the segmentation function 
#   @seemethod "DNAcopy::segment" is used.
#   @seeclass
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("segmentByCBS", "RawGenomicSignals", function(this, ..., verbose=FALSE) {
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
  pkgName <- "DNAcopy";
  # Assert that package is installed
  isPackageInstalled(pkgName) || throw("Package is not installed: ", pkgName);
  pkg <- packageDescription(pkgName);
  pkgVer <- pkg$Version;
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer);

  methodName <- "segment";
  verbose && cat(verbose, "Method: ", methodName);
  verbose && cat(verbose, "Package: ", pkgDetails);

  # We need to load package
  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);

  # Get the fit function for the segmentation method
  fitFcn <- getExportedValue(pkgName, methodName);
  fitFcn <- getFromNamespace(methodName, pkgName);
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
  hasWeights <- !is.null(data$weight);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasWeights) {
    # Verify that weights are supported (from DNAcopy v1.18.0);
    if (!is.element("weights", names(formals))) {
      hasWeights <- FALSE;
      msg <- paste("Weights detected but ignored, because the available segmentation function does not support weights. Check with a more recent version of the package: ", pkgDetails);
      verbose && cat(verbose, msg);
      warning(msg);
    }
  } # if (hasWeights)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure");
  cnData <- DNAcopy::CNA(
    genomdat  = data$signal, 
    chrom     = data$chromosome,
    data.type = "logratio",
    maploc    = data$position,
    sampleid  = sampleName
  );
  verbose && str(verbose, cnData);
  verbose && exit(verbose);

  args <- list(cnData);

  if (hasWeights) {
    fitArgs <- list(weights=data$weight);
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
  # Note: DNAcopy::segment() *does* this.
  stdout <- capture.output({
    # Does not work, because some internal function of the fit function
    # may only be accessible from within the namespace
    # How to do this for DNAcopy::segment()? /HB
##    fit <- do.call(fitFcn, args);
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    fit <- do.call(methodName, args);
  });

  verbose && cat(verbose, "Captured output that was sent to stdout:");
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

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
