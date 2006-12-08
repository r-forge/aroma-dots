setMethodS3("importFeatureSet", "AffymetrixCelFile", function(this, featureSet, ..., .map=NULL, unprotect=FALSE, .checkArgs=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    if (!is(featureSet, "FeatureSet")) {
      throw("Argument 'featureSet' is not a FeatureSet object: ", 
                                                           class(featureSet)[1]);
    }
    dim <- dim(featureSet);
    if (dim[2] < 0) {
      throw("Argument 'featureSet' contains no data.");
    } else if (dim[2] > 1) {
      throw("Argument 'featureSet' must only contain on sample: ", dim[2]);
    }

    if (!identical(unprotect, TRUE)) {
      throw("In order to import data to a CEL file, set argument 'unprotect' to TRUE.");
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing probe signals to CEL file");
  pathname <- getPathname(this);
  verbose && cat(verbose, "Sample name: ", getName(this));
  verbose && cat(verbose, "Tags: ", paste(getTags(this), collapse=","));
  verbose && cat(verbose, "Pathname: ", pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate the feature set  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    verbose && enter(verbose, "Validating that the import data is compatible with the chip type");
    cdf <- getCdf(this);
    cdf2 <- AffymetrixCdfFile$fromCleanChipType(annotation(featureSet));
    if (!identical(getChipType(cdf), getChipType(cdf2))) {
      throw("The CDF of the feature set data is not compatible with the CEL file: ", getChipType(cdf), " != ", getChipType(cdf2));
    }
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the platform-design feature index to CDF cell index map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(.map)) {
    cdf <- getCdf(this);
    .map <- getOligoToCelMap(cdf, verbose=less(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Storing probe signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Storing probe intensities");
  y <- unname(assayData(featureSet)$exprs[.map,1]);
  updateCel(pathname, intensities=y, verbose=-as.integer(verbose)-5);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(.map);  
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-12-06
# o Created.
############################################################################
