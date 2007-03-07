setMethodS3("setSampleAttributesByTags", "AffymetrixCelSet", function(this, sampleName, tags, ...) {
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);

  # Identify all files that refers to this sample
  idxs <- which(sampleName == getNames(this));
  if (length(idxs) == 0)
    return();

  # Set the attribute for all such files
  newAttrs <- list();
  for (idx in idxs) {
    cf <- getFile(cs, idx);
    newAttrs <- c(newAttrs, setAttributesByTags(cf, tags, ...));
  }

  invisible(newAttrs);
}, protected=TRUE)

############################################################################
# HISTORY:
# 2007-03-06
# o Added setSampleAttributesByTags().
# o Created.
############################################################################
