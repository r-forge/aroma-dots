setConstructorS3("SdrfFileSet", function(...) {
  extend(TabularTextFileSet(...), "SdrfFileSet");
})


setMethodS3("byPath", "SdrfFileSet", function(static, ..., pattern=".*[.]sdrf[.]txt$") {
  suppressWarnings({
    byPath.GenericDataFileSet(static, ..., pattern=pattern);
  })
})


## setMethodS3("fromFileSet", "GenericDataFileSet", function(static, set, ...) {
##   # Argument 'set':
##   className <- "GenericDataFileSet"
##   if (!inherits(set, className)) {
##     throw("Argument 'set' is not a ", className, ": ", class(set)[1]);
##   }
##  
##   path <- getPath(set);
##   fromFiles(static, path=path, ...);
## })



setMethodS3("getFilenameMap", "SdrfFileSet", function(this, ...) {
  # Extract the data for each file
  data <- lapply(this, FUN=function(file) {
    data <- getFilenameMap(file, ...);
    # Drop the rownames temporarily
    data$fullName <- rownames(data);
    rownames(data) <- NULL;
    data;
  });
  
  # Stack them
  data <- Reduce(rbind, data);

  # Drop duplicates
  keep <- which(!duplicated(data$arrayDataFile));
  data <- data[keep,,drop=FALSE];

  # Add rownames
  rownames(data) <- data$fullName;

  data$fullName <- NULL;

  data;
})


setMethodS3("makeFullNamesTranslator", "SdrfFileSet", function(this, ..., failTag="notInSDRF", force=FALSE) {
  translator <- this$.fullnamesTranslator;

  if (force || is.null(translator)) {
    # Extract the data for each file
    data <- getFilenameMap(this, ...);
  
    translator <- function(names, ...) {
      names2 <- data[names,"extractName"];
      # Any failed translations?
      nok <- is.na(names2);
      names2[nok] <- paste(names[nok], failTag, sep=",");
      names2;
    } # translator()

    this$.fullnamesTranslator <- translator;
  }

  # Return translator function
  translator;
})



############################################################################
# HISTORY:
# 2009-10-30
# o Now makeFullNamesTranslator() of SdrfFile|Set adds a "failTag" to the
#   fullname of those fullnames that was not translated.
# o Added getFilenameMap() to SdrfFileSet.
# 2009-10-23
# o OPTIMIZATION: Now makeFullNamesTranslator() of SdrfFileSet caches the
#   result in the memory so it does not have to reread the files each times.
# 2009-10-02
# o Created.
############################################################################
