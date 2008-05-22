###########################################################################/**
# @RdocClass ScpModel
#
# @title "The ScpModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Stochastic Change-Point (SCP) model [1]. 
#  The required package is currently private, but according to the
#  authors, it will soon be released.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
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
#
# \references{
#  [1] Lai et al., \emph{Stochastic segmentation models for array-based 
#      comparative genomic hybridization data analysis}, 
#      Biostatistics, 2008.\cr
# }
#*/###########################################################################
setConstructorS3("ScpModel", function(cesTuple=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    # Fool R CMD check (because 'cnv' is not publically available)
    pkgName <- "cnv";
    require(pkgName, character.only=TRUE) || throw("Package not loaded: cnv");
  }

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "ScpModel")
})


setMethodS3("getAsteriskTags", "ScpModel", function(this, collapse=NULL, ...) {
  tags <- "SCP";

  # Add class-specific tags
  if (isPaired(this))
    tags <- c(tags, "paired");

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("fitOne", "ScpModel", function(this, data, chromosome, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Fitting SCP");

  verbose && enter(verbose, "Estimating raw CN variance");
  s2 <- varDiff(data[,"M"], na.rm=TRUE);
  verbose && printf(verbose, "s2 = %.2g", s2);
  verbose && exit(verbose);

  # Default hyper parameters (from example script)
  hyperParams <- list(p=0.95, a=0.99, b=0.005, mu=0, v=1, sigma2=s2);

  # Parameters (from example script)
  params <- list(hyper=hyperParams, classify="T", thres=0.1);
  verbose && cat(verbose, "Default paramters:");
  verbose && str(verbose, params);

  # Override with user specified parameters
  args <- list(...);
  # Fool R CMD check (because we do not wanna use cnv::... yet)
  getBcmixSmoothClass <- NULL; rm(getBcmixSmoothClass);
  keep <- (names(args) %in% names(formals(getBcmixSmoothClass)));
  args <- args[keep];
  for (key in names(args)) {
    params[[key]] <- args[[key]];
  }
  verbose && cat(verbose, "Paramters overridden by user arguments:");
  verbose && str(verbose, params);


  verbose && enter(verbose, "Setting up SCP data structure");
  nbrOfUnits <- nrow(data);
  chipTypes <- getChipTypes(this);

  xPos <- data[,"x"];
  data <- list(
    obs = data[,"M"]
  );
  verbose && str(verbose, data);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling getBcmixSmoothClass()");
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits);
  args <- c(data, params);
  verbose && cat(verbose, "All arguments to fit function:");
  verbose && str(verbose, args);
  rm(data, params);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # getBcmixSmoothClass() writes to stdout; capture it and send it
  # to the verbose object.
  stdout <- capture.output({
    fit <- do.call("getBcmixSmoothClass", args);
  })
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);
  verbose && str(verbose, fit);

  verbose && exit(verbose);

  verbose && enter(verbose, "Adding physical position to the fit object");
  fit$xPos <- xPos;
  fit$chromosome <- chromosome;
  class(fit) <- c("SCPfit", class(fit));
  verbose && str(verbose, fit);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;  
}, private=TRUE) # fitOne()


setMethodS3("extractRawCopyNumbers", "SCPfit", function(object, ...) {
  RawCopyNumbers(cn=object$obs, x=object$xPos, chromosome=object$chromosome);
})


setMethodS3("extractCopyNumberRegions", "SCPfit", function(object, ...) {
  CopyNumberRegions();
})



##############################################################################
# HISTORY:
# 2008-05-21
# o Now extractRawCopyNumbers() adds 'chromosome' to the returned object.
# 2008-04-17.
# o Created.
##############################################################################
