# path <- file.path("rawData", getName(ds));
# sa <- SampleAnnotationFile("sampleAnnotation.xls", path=path);

setConstructorS3("SampleAnnotationFile", function(...) {
  this <- extend(AffymetrixFile(...), "SampleAnnotationFile",
    "cached:.data" = NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    parseTagsAsAttributes(this);

  this;
})



setMethodS3("getAlias", "SampleAnnotationFile", function(this, ...) {
  # Get the records of interest.
  data <- getDataByName(this, ...);
  names <- rownames(data);
  data <- data[,"alias"];

  # Non-listed samples
  nas <- which(is.na(data));
  data[nas] <- names[nas];

  data;
})



setMethodS3("getNbrOfChrX", "SampleAnnotationFile", function(this, ...) {
  # Get the records of interest.
  data <- getDataByName(this, ...);
  data <- data[,"nbrOfChrX"];

  # Convert to numeric
  data <- as.double(data);

  data;
})



setMethodS3("getDataByName", "SampleAnnotationFile", function(this, names=NULL, ...) {
  # Argument 'names':
  names <- Arguments$getCharacters(names);

  # Get all annotation data
  data <- getData(this);

  if (is.null(names))
    return(data);

  # Locate the rows for each of the requested samples
  patterns <- paste("^", rownames(data), "$", sep="");
  matches <- sapply(patterns, FUN=grep, names, value=TRUE);

  # Find each of the samples
  rr <- integer(length(names));
  for (kk in seq(along=names)) {
    name <- names[kk];
    idxs <- lapply(matches, FUN=function(xs) { name %in% xs });
    idxs <- which(unlist(idxs, use.names=FALSE));
    rr[kk] <- idxs[1];
  }

  data <- data[rr,,drop=FALSE];
  rownames(data) <- names;

  data;
}, private=TRUE)



setMethodS3("getData", "SampleAnnotationFile", function(this, ...) {
  data <- this$.data;
  if (is.null(data)) {
    data <- readData(this);
    this$.data <- data;
  }
  data;
})



setMethodS3("readData", "SampleAnnotationFile", function(this, ...) {
  pathname <- getPathname(this);
  df <- readTable(pathname, header=TRUE, na.strings=c("NA", "", "-"), ...);

  cc <- match("name", colnames(df));
  rownames(df) <- df[,cc];
  df <- df[,-cc,drop=FALSE];

  df;
}, private=TRUE)


############################################################################
# HISTORY:
# 2007-01-26
# o Created.
############################################################################
