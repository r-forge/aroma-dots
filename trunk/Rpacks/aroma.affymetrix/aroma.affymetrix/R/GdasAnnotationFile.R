setConstructorS3("GdasAnnotationFile", function(...) {
  extend(AffymetrixFile(...), "GdasAnnotationFile",
    "cached:.blockPositions" = NULL,
    "cached:.header" = NULL,
    "cached:.columnName" = NULL,
    "cached:.data" = NULL
  )
})


setMethodS3("clearCache", "GdasAnnotationFile", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".blockPositions", ".header", ".columnName", ".data")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)



setMethodS3("fromChipType", "GdasAnnotationFile", function(static, chipType, what, path="annotations", ...) {
  # Argument 'chipType':
  if (inherits(chipType, "AffymetrixCdfFile")) {
    chipType <- getChipType(chipType);
  } else if (is.character(chipType)) {
  } else {
    throw("Argument 'chipType' must be a character string or an AffymetrixCdfFile object: ", class(chipType)[1]);
  }

  fullname <- sprintf("%s_%s", chipType, what);
  filename <- sprintf("%s.tsv", fullname);

  newInstance(static, filename=filename, path=path, ...);
}, static=TRUE)


setMethodS3("getBlockPositions", "GdasAnnotationFile", function(this, ..., force=FALSE) {
  # Check cache
  blockPositions <- this$.blockPositions;
  if (!force && !is.null(blockPositions))
    return(blockPositions);

  pathname <- getPathname(this);
  lines <- readLines(pathname);
  nlines <- length(lines);
  pattern <- "^\\[(.*)\\]";
  pos <- grep(pattern, lines);
  lines <- lines[pos];
  names <- gsub(pattern, "\\1", lines);
  pos <- matrix(c(pos+1,pos[-1]-1,nlines), ncol=2);
  colnames(pos) <- c("from", "to");
  rownames(pos) <- names;

  # Save to cache
  this$.blockPositions <- pos;

  pos;
}, private=TRUE)


setMethodS3("readBlock", "GdasAnnotationFile", function(this, name, ...) {
  # Find the position of the block
  pos <- getBlockPositions(this, ...);
  if (!name %in% rownames(pos))
    throw("Cannot read block. No such block: ", name);
  pos <- pos[name,];

  # Read the block
  pathname <- getPathname(this);
  con <- file(pathname, open="r");
  on.exit(close(con));

  # Skip to block
  readLines(con, n=pos-1);

  # Read block
  lines <- readLines(con, n=diff(pos)+1);

  lines;
}, private=TRUE)


setMethodS3("readHeader", "GdasAnnotationFile", function(this, ...) {
  header <- this$.header;
  if (!is.null(header))
    return(header);

  # Read from file
  res <- readBlock(this, "Header");

  # Parse it
  names <- c();
  header <- vector("list", length(res));
  for (kk in seq(along=res)) {
    value <- res[kk];
    value <- trim(value);
    pos <- regexpr("=", value)[1];
    name <- substring(value, 1, pos-1);
    name <- trim(name);
    value <- substring(value, pos+1);
    value <- trim(value);
    names[kk] <- name;
    header[[kk]] <- value;
  }
  names(header) <- names;

  # Save to cached field
  this$.header <- header;

  header;
})


setMethodS3("readColumnName", "GdasAnnotationFile", function(this, ..., force=FALSE) {
  columnName <- this$.columnName;
  if (!force && !is.null(columnName))
    return(columnName);

  # Read from file
  res <- readBlock(this, "ColumnName");

  # Parse it
  columnName <- strsplit(res, split="\t")[[1]];

  # Save to cached field
  this$.columnName <- columnName;

  columnName;
})



setMethodS3("readData", "GdasAnnotationFile", function(this, colClasses="character", ..., force=FALSE) {
  data <- this$.data;
  if (!force && !is.null(data))
    return(data);

  # Read from file
  data <- readBlock(this, "Data");
  data <- strsplit(data, split="\t");

  colnames <- readColumnName(this, force=force);
  ncol <- length(colnames);

  data <- unlist(data, use.names=FALSE);
  data <- matrix(data, nrow=ncol);
  data <- t(data);
  colnames(data) <- colnames;
  data <- as.data.frame(data, stringsAsFactors=FALSE);

  colClasses <- rep(colClasses, length.out=ncol);
  for (cc in seq(length=ncol)) {
    value <- data[,cc];
    storage.mode(value) <- colClasses[cc];
    data[,cc] <- value;
  }

  # Save to cached field
  this$.data <- data;

  data;
})


setMethodS3("[", "GdasAnnotationFile", function(this, i=NULL, j=NULL, drop=TRUE) {
  data <- readData(this);
  if (is.character(i)) {
    i <- match(i, data[,1]);
  }
  if (!is.null(i))
    data <- data[i,,drop=FALSE];
  if (!is.null(j))
    data <- data[,j,drop=FALSE];
  if (drop)
    data <- drop(data);
  data;
})


############################################################################
# HISTORY:
# 2006-08-27
# o Created.
############################################################################
