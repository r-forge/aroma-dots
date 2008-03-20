setMethodS3("getAsFullCelSet", "ChipEffectSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  verbose && enter(verbose, "Getting chip effect set expanded to the full CDF");
  files <- list();
  for (kk in seq(this)) {
    cef <- getFile(this, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') %d", 
                                          kk, getName(cef), length(this)));
    cf <- getAsFullCelFile(cef, ..., verbose=less(verbose, 5));
    files[[kk]] <- cf;
  }
  verbose && exit(verbose);

  res <- AffymetrixCelSet(files);

  res;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2008-03-20
# o Created.
############################################################################
