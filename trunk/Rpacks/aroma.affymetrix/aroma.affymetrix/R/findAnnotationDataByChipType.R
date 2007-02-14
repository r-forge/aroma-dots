setMethodS3("findAnnotationDataByChipType", "default", function(chipType, pattern=chipType, ..., paths=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get search paths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get paths to search
  if (is.null(paths)) {
    settings <- getOption("aroma.affymetrix.settings");
    paths <- settings$paths$annotationData;
  }

  if (is.null(paths)) {
    paths <- "annotationData";
  } else {
    # Split path strings by semicolons.
    paths <- unlist(strsplit(paths, split=";"));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Expand any file system links
  paths <- file.path(paths, "chipTypes", chipType);
  paths <- sapply(paths, FUN=filePath, expandLinks="any");

  # Search recursively for all CDF files
  pathname <- findFiles(pattern, paths=paths, recursive=TRUE, ...);

  pathname;
}, protected=TRUE)

############################################################################
# HISTORY:
# 2007-02-06
# o Created.
############################################################################
