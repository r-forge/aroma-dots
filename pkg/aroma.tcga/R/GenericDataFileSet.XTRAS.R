setMethodS3("appendFullNameTranslatorBySdrfFile", "FullNameInterface", function(this, sdrf, ...) {
  sdrfs <- SdrfFileSet(sdrf);
  appendFullNameTranslator(this, sdrfs);
})


setMethodS3("appendFullNameTranslatorBySdrfFileSet", "FullNameInterface", function(this, set, ...) {
  # Argument 'set':
  set <- Arguments$getInstanceOf(set, "SdrfFileSet");

  # Build a fullnames translator function based on the SDRF
  fcn <- makeFullNamesTranslator(set, ...);

  # Apply the translator function
  appendFullNameTranslator(this, fcn);
})


setMethodS3("extractByTcgaType", "GenericDataFileSet", function(this, pattern, ..., maxNbrOfNonTcgaSamples=0.1) {
  getTypeTag <- function(this, ...) {
    gsub("[A-Z]$", "", getTags(this)[1]);
  }

  # Argument 'typePattern':
  pattern <- Arguments$getRegularExpression(pattern);

  # Argument 'maxNbrOfNonTcgaSamples':

  # Sanity check
  names <- getFullNames(this);

  if (maxNbrOfNonTcgaSamples > 0 && maxNbrOfNonTcgaSamples < 1) {
    maxNbrOfNonTcgaSamples <- round(maxNbrOfNonTcgaSamples*length(names));
  }

  patterns <- BiospecimenCoreResource$getBarcodePatterns();
  pattern2 <- sprintf("^%s,%s", patterns$patient, patterns$sampleId);
  bad <- which(regexpr(pattern2, names) == -1);
  if (length(bad) > maxNbrOfNonTcgaSamples) {
    throw(sprintf("%s '%s' contains %s files with non-TCGA names: %s", 
          class(this)[1], getFullName(this), length(bad)), 
          paste(bad, collapse=", "));
  }

  # Identify sample types matching the pattern
  types <- sapply(this, getTypeTag);
  keep <- grep(pattern, types);

  # Extract matching samples
  res <- extract(this, keep);

  res;
})


setMethodS3("extractTumorNormalPairs", "GenericDataFileSet", function(this, tumorPattern="^(01|20)", normalPattern="^(10|11)", ..., onDuplicates=c("drop", "ignore"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tumorPattern':
  tumorPattern <- Arguments$getRegularExpression(tumorPattern);

  # Argument 'normalPattern':
  normalPattern <- Arguments$getRegularExpression(normalPattern);

  # Argument 'onDuplicates':
  onDuplicates <- match.arg(onDuplicates);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting (tumor, normal) pairs");
  verbose && cat(verbose, "Tumor type tag pattern: ", tumorPattern);
  verbose && cat(verbose, "Normal type tag pattern: ", normalPattern);

  verbose && enter(verbose, "Identifying tumors and normals");
  dsT <- extractByTcgaType(this, pattern=tumorPattern);
  dsN <- extractByTcgaType(this, pattern=normalPattern);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying common set");
  # Keep only samples for which there is both a tumor and a normal
  sampleNames <- intersect(getNames(dsT), getNames(dsN));
  sampleNames <- unique(sampleNames);
  verbose && cat(verbose, "Identified tumor-normal samples:");
  verbose && print(verbose, sampleNames);
  verbose && exit(verbose);

  dsT <- extract(dsT, sampleNames);
  dsN <- extract(dsN, sampleNames);

  # Drop duplicated samples
  if (onDuplicates == "drop") {
    verbose && enter(verbose, "Dropping duplicates");
    dsT <- extract(dsT, !duplicated(getNames(dsT)));
    dsN <- extract(dsN, !duplicated(getNames(dsN)));
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Data set of tumors:");
  verbose && print(verbose, dsT);
  verbose && cat(verbose, "Data set of normals:");
  verbose && print(verbose, dsN);

  verbose && exit(verbose);

  list(tumor=dsT, normal=dsN);
})



############################################################################
# HISTORY:
# 2010-01-08
# o Added extractTumorNormalPairs().
# 2010-01-05 
# o Added extractByTcgaType().
# 2009-10-23
# o OPTIMIZATION: Now makeFullNamesTranslator() of SdrfFileSet caches the
#   result in the memory so it does not have to reread the files each times.
# 2009-10-02
# o Added setFullNamesTranslatorBySdrfFileSet() to GenericDataFileSet.
# 2009-10-01
# o Added setFullNamesTranslatorBySdrfFile() to GenericDataFileSet.
# o Created.
############################################################################
