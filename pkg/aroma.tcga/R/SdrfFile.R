setConstructorS3("SdrfFile", function(..., .verify=TRUE) {
  this <- extend(TabularTextFile(..., .verify=FALSE), "SdrfFile");

  if (!is.null(getPathname(this))) {
    setColumnNameTranslator(this, toCamelCase);
  }
#  if (.verify)
#    verify(this, ..., verbose=verbose);
  this;
})


setMethodS3("getFilenameMap", "SdrfFile", function(this, ...) {
  colClassPatterns <- c("(extractName|arrayDataFile)"="character");
  readDataFrame(this, colClassPatterns=colClassPatterns);
})


setMethodS3("getExtractNameByFilename", "SdrfFile", function(this, name, ...) {
  colClassPatterns <- c("(extractName|arrayDataFile)"="character");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns);
  arrayDataFile <- NULL; rm(arrayDataFile); # To trick R CMD check
  data <- subset(data, arrayDataFile == name);
  data$extractName;
})


setMethodS3("makeFullNamesTranslator", "SdrfFile", function(this, ...) {
  data <- getFilenameMap(this, ...);
  data$fullNames <- gsub("[.](CEL|cel)$", "", data$arrayDataFile);

  # TO DO: The translator function must accept test names 'foo' and 'bar'.
  safe <- FALSE;

  translator <- function(names, ...) {
    idxs <- match(names, data$fullNames);

    # Are all fullnames defined?
    ok <- is.finite(idxs);
    if (safe && !all(ok)) {
      missing <- names[!ok];
      missing <- paste(head(missing), collapse=", ");
      throw("Failed to translate the name of some data files. Unknown fullnames: ", missing);
    }

    names[ok] <- data$extractName[idxs[ok]];
    names;
  } # translator()

  # Return translator function
  translator;
})


############################################################################
# HISTORY:
# 2009-10-01
# o Added makeFullNamesTranslator() to SdrfFile.
# 2009-05-21
# o Created.
############################################################################
