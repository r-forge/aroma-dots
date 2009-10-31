setConstructorS3("SdrfFile", function(..., .verify=TRUE) {
  this <- extend(TabularTextFile(..., .verify=FALSE), "SdrfFile");

  if (!is.null(getPathname(this))) {
    setColumnNameTranslator(this, toCamelCase);
  }
#  if (.verify)
#    verify(this, ..., verbose=verbose);
  this;
})


setMethodS3("getFilenameMap", "SdrfFile", function(this, ..., dropEmpty=TRUE, .validate=FALSE) {
  colClassPatterns <- c("(extractName|arrayDataFile)"="character");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns);

  # Trim (just in case)
  for (key in colnames(data)) {
    data[,key] <- trim(data[,key]);
  }    

  # Drop duplicates
  keep <- which(!duplicated(data$arrayDataFile));
  data <- data[keep,,drop=FALSE];

  for (key in colnames(data)) {
    data[,key] <- trim(data[,key]);
  }    

  # Validate
  if (.validate) {
    for (key in colnames(data)) {
      values <- data[,key];
      len <- nchar(values);
      idxs <- which(len == 0);
      if (length(idxs) > 0) {
        msg <- sprintf("Invalid filename map. Column '%s' contains empty values: Rows %s", key, paste(head(idxs), collapse=", "));
        throw(msg);
      }
    } # for (key ...)  
  }

  # Drop empty target names?
  if (dropEmpty) {
    keep <- which(nchar(data$extractName) > 0);
    data <- data[keep,,drop=FALSE];
  }

  # Add rownames
  # Drop the filename extension
  pattern <- "\\.[a-zA-Z]+$";
  pattern <- toAsciiRegExprPattern(pattern);
  rownames(data) <- gsub(pattern, "", data$arrayDataFile);

  data;
})


setMethodS3("getExtractNameByFilename", "SdrfFile", function(this, name, ...) {
  colClassPatterns <- c("(extractName|arrayDataFile)"="character");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns);
  arrayDataFile <- NULL; rm(arrayDataFile); # To trick R CMD check
  data <- subset(data, arrayDataFile == name);
  data$extractName;
})


setMethodS3("makeFullNamesTranslator", "SdrfFile", function(this, ..., failTag="notInSDRF") {
  data <- getFilenameMap(this, ...);

  translator <- function(names, ...) {
    names2 <- data[names,"extractName"];
    # Any failed translations?
    nok <- is.na(names2);
    names2[nok] <- paste(names[nok], failTag, sep=",");
    names2;
  } # translator()

  # Return translator function
  translator;
})


############################################################################
# HISTORY:
# 2009-10-30
# o Now makeFullNamesTranslator() of SdrfFile|Set adds a "failTag" to the
#   fullname of those fullnames that was not translated.
# o Now the data frame returned by getFilenameMap() of SdrfFile has the
#   default fullnames as rownames.  It also drops all rows that have empty
#   extractName:s.
# o ROBUSTNESS: Now getFilenameMap() of SdrfFile can assert that there are
#   no empty values.
# 2009-10-23
# o Now makeFullNamesTranslator() of SdrfFile translates all data types.
# 2009-10-01
# o Added makeFullNamesTranslator() to SdrfFile.
# 2009-05-21
# o Created.
############################################################################
