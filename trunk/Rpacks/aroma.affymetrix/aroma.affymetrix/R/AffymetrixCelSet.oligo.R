setMethodS3("createFromFeatureSet", "AffymetrixCelSet", function(static, featureSet, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is(featureSet, "FeatureSet")) {
    throw("Argument 'featureSet' is not a FeatureSet object: ", 
                                                         class(featureSet)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the sample names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sampleNames <- sampleNames(phenoData(featureSet));
  # Remove any .CEL at the end
  sampleNames <- gsub("[.](c|C)(e|E)(l|L)$", "", sampleNames);
  nbrOfSamples <- length(sampleNames);

  verbose && enter(verbose, "Creating ", class(static)[1], " from ", class(featureSet)[1]);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create set of blank CEL files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- createBlankSet(static, ..., sampleNames=sampleNames, verbose=verbose);

  map <- NULL;
  for (kk in seq(ds)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array #%d of %d", kk, nbrOfSamples));
    map <- importFeatureSet(df, featureSet=featureSet[,kk], .checkArgs=(kk==1), .map=map, unprotect=TRUE, verbose=verbose);
    verbose && exit(verbose);
  }
  verbose && exit(verbose);
  
  ds;
}, static=TRUE);




############################################################################
# HISTORY:
# 2006-12-06
# o Created.
############################################################################
