setConstructorS3("TextUnitNamesFile", function(..., platform=NULL) {
  # Argument 'platform':
  if (!is.null(platform)) {
    platform <- Arguments$getCharacter(platform);
  }

  extend(TabularTextFile(...), c("TextUnitNamesFile", uses("UnitNamesFile")),
    .platform = platform
  );
})


setMethodS3("as.character", "TextUnitNamesFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", this, ...);
  class <- class(s);
  s <- c(s, sprintf("Number units: %d", nbrOfUnits(this, fast=TRUE)));

  class(s) <- class;
  s;
})

setMethodS3("getPlatform", "TextUnitNamesFile", function(this, ...) {
  platform <- this$.platform;
  if (is.null(platform)) {
    comments <- getHeader(this)$comments;
    params <- gsub("^#", "", comments);
    params <- trim(params);
    params <- strsplit(params, split=":");
    params <- lapply(params, FUN=trim);
    keys <- sapply(params, FUN=function(x) x[1]);
    values <- sapply(params, FUN=function(x) x[-1]);
    names(values) <- keys;
    platform <- unname(values["platform"]);
    if (is.na(platform))
      platform <- NULL;
    this$.platform <- platform;
  }
  platform;
})

setMethodS3("getChipType", "TextUnitNamesFile", function(this, ...) {
  getName(this);
})

setMethodS3("getUnitNames", "TextUnitNamesFile", function(this, units=NULL, ...) {
  data <- readDataFrame(this, rows=units, ...);
  unitNames <- data[,1,drop=TRUE];
  unitNames <- as.character(unitNames);
  if (is.null(units)) {
    this$.nbrOfUnits <- length(unitNames);
  }
  unitNames;
})

setMethodS3("nbrOfUnits", "TextUnitNamesFile", function(this, ...) {
  nbrOfUnits <- this$.nbrOfUnits;
  if (is.null(nbrOfUnits)) {
    nbrOfUnits <- NextMethod("nbrOfUnits", this);
    this$.nbrOfUnits <- nbrOfUnits;
  }
  nbrOfUnits;
})

setMethodS3("getFilenameExtension", "TextUnitNamesFile", function(static, ...) {
  "txt";
}, static=TRUE, protected=TRUE);

setMethodS3("findByChipType", "TextUnitNamesFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get fullname, name, and tags
  fullname <- paste(c(chipType, tags), collapse=",");
  parts <- unlist(strsplit(fullname, split=","));
  # Strip 'monocell' parts
  parts <- parts[parts != "monocell"];
  chipType <- parts[1];
  tags <- parts[-1];
  fullname <- paste(c(chipType, tags), collapse=",");

  ext <- getFilenameExtension(static);
  ext <- paste(c(tolower(ext), toupper(ext)), collapse="|");
  ext <- sprintf("(%s)", ext);

  pattern <- sprintf("^%s.*,unitNames[.]%s$", fullname, ext);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;  # Override argument 'pattern'?
#  args$firstOnly <- FALSE;
#  str(args);
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- sprintf("^%s.*[.]%s[.]lnk$", chipType, ext);
    args$pattern <- pattern;
    pathname <- do.call("findAnnotationDataByChipType", args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- filePath(pathname, expandLinks="any");
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2009-04-06
# o BUG FIX: getUnitNames(..., units=NULL) of TextUnitNamesFile would 
#   make the object believe there are zero units in the file.
# 2009-03-23
# o Created.
############################################################################
