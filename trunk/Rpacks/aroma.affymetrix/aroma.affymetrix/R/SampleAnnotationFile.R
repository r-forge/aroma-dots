# path <- getParent(getPath(cs))
# sa <- SampleAnnotationFile("sampleAnnotation.dcf", path=path);

setConstructorS3("SampleAnnotationFile", function(...) {
  this <- extend(AffymetrixFile(...), "SampleAnnotationFile",
    "cached:.db" = NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})

setMethodS3("readData", "SampleAnnotationFile", function(this, rows=NULL, force=FALSE, ...) {
  db <- this$.db;
  if (force || is.null(db)) {
    pathname <- getPathname(this);
  
    # Read all non-commented lines
    bfr <- readLines(pathname); 
    bfr <- bfr[-grep("^[ ]*#", bfr)];
  
    # Parse these as a DCF
    db <- read.dcf(textConnection(bfr));
    rm(bfr); # Not needed anymore
  
    this$.db <- db;
  }

  colnames(db) <- toCamelCase(colnames(db));

  if (!is.null(rows))
    db <- db[rows,,drop=FALSE];

  db;
}, protected=TRUE)


setMethodS3("getPatterns", "SampleAnnotationFile", function(this, ...) {
  db <- readData(this, ...);

  # Get sample name pattern
  patterns <- sprintf("^%s$", db[,"sample"]);
  patterns <- gsub("\\^\\^", "^", patterns);
  patterns <- gsub("\\$\\$", "$", patterns);

  patterns;
}, protected=TRUE)

setMethodS3("match", "SampleAnnotationFile", function(this, names, trim=FALSE, ...) {
  # Scan vector of names for matching patterns
  patterns <- getPatterns(this, ...);
  res <- sapply(patterns, FUN=function(pattern) { 
    idxs <- grep(pattern, names);
    names(idxs) <- names[idxs];
    idxs;
  });

  if (trim) {
    keep <- (sapply(res, FUN=length) > 0);
    res <- res[keep];
  }

  res;
}, protected=TRUE)


setMethodS3("apply", "SampleAnnotationFile", function(this, names, FUN, ...) {
  allPatterns <- getPatterns(this, ...);
  res <- match(this, names, trim=TRUE);
  patterns <- names(res);
  rows <- match(patterns, allPatterns);
  db <- readData(this, rows=rows);
  cc <- setdiff(colnames(db), "sample");
  db <- db[,cc,drop=FALSE];

  for (kk in seq(along=res)) {
    record <- db[kk,,drop=TRUE];

    # Nothing to do?
    if (all(is.na(record)))
      next;

    args <- list(
      appliesTo = res[[kk]]
    );
    args <- c(args, as.list(record));
    args <- c(args, list(...));
    do.call("FUN", args=args);
  }
})

############################################################################
# HISTORY:
# 2007-01-26
# o Created.
############################################################################
