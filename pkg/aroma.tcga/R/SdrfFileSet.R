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



setMethodS3("makeFullNamesTranslator", "SdrfFileSet", function(this, ..., force=FALSE) {
  translator <- this$.fullnamesTranslator;

  if (force || is.null(translator)) {
    # Extract the data for each file
    data <- lapply(this, getFilenameMap, ...);
  
    # Stack them
    data <- Reduce(rbind, data);
  
    # Translate
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

    this$.fullnamesTranslator <- translator;
  }

  # Return translator function
  translator;
})



############################################################################
# HISTORY:
# 2009-10-23
# o OPTIMIZATION: Now makeFullNamesTranslator() of SdrfFileSet caches the
#   result in the memory so it does not have to reread the files each times.
# 2009-10-02
# o Created.
############################################################################
