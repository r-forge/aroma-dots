setConstructorS3("GdasAnnotationSet", function(...) {
  extend(AffymetrixFileSet(...), "GdasAnnotationSet",
    .cdf = NULL
  )
})

setMethodS3("fromFiles", "GdasAnnotationSet", function(static, path="annotations", ...) {
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  pattern <- "[.]tsv$";
  pathnames <- findFiles(pattern=pattern, path=path, firstOnly=FALSE);
  if (length(pathnames) == 0)
    throw("No files found: ", path);

  files <- lapply(pathnames, FUN=GdasAnnotationFile);

  newInstance(static, files);
}, static=TRUE)

setMethodS3("fromChipType", "GdasAnnotationSet", function(static, chipType, what, ...) {
  files <- vector("list", length(what));
  for (kk in seq(along=what)) {
    files[[kk]] <- GdasAnnotationFile$fromChipType(chipType=chipType, what=what[kk], ...);
  }

  newInstance(static, files);
}, static=TRUE)

setMethodS3("join", "GdasAnnotationSet", function(this, indexColumn=1, orderBy=NULL, ..., force=FALSE) {
  df <- this$.join;
  if (force || is.null(df)) {
    df <- NULL;
    files <- as.list(this);
    for (file in files) {
      data <- readData(file);
      if (is.null(df)) {
        df <- data;
      } else {
        df <- merge(df, data);
      }
    }
    this$.join <- df;
  }

  if (!is.null(orderBy)) {
    args <- vector("list", length(orderBy));
    for (kk in seq(along=args))
      args[[kk]] <- df[,orderBy[kk]];
    o <- do.call("order", args=args);
    df <- df[o,];
    rownames(df) <- 1:nrow(df);
  }
  this$.join <- df;

  df;
})


setMethodS3("[", "GdasAnnotationSet", function(this, i=NULL, j=NULL, drop=TRUE) {
  data <- join(this);
  names <- data[,1];
  data <- data[,-1];
  rownames(data) <- names;
  if (!is.null(i))
    data <- data[i,,drop=FALSE];
  if (!is.null(j))
    data <- data[,j,drop=FALSE];
  data;
})


setMethodS3("getChipType", "GdasAnnotationSet", function(this, ...) {
  # Figure out the chip type from the 'ArrayName' header in the files
  file <- getFile(this, 1);
  hdr <- readHeader(file);
  chipType <- hdr$ArrayName;

  # Check if we can find a CDF for this chip type
  cdfFile <- findCdf(chipType);
  if (is.null(cdfFile))
    warning("Could not find a CDF file for inferred chip type: ", chipType);
    
  chipType;
})

setMethodS3("getCdf", "GdasAnnotationSet", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    cdf <- AffymetrixCdfFile$fromChipType(getChipType(this));
    this$.cdf <- cdf;
  }
  cdf;
})

setMethodS3("select", "GdasAnnotationSet", function(this, ..., sortBy=NULL) {
  dfAll <- this[];

  args <- list(...);
  fields <- names(args);

  df <- dfAll;
  for (kk in seq(along=fields)) {
    field <- fields[kk];
    field <- strsplit(field, split=":")[[1]];
    if (length(field) > 1) {
      modifier <- field[1];
      field <- paste(field[-1], collapse=":");
    } else {
      modifier <- NULL;
    }

    # Partial matching for field names
    cc <- pmatch(field, colnames(df));

    data <- df[,cc];
    test <- args[[kk]];
    if (is.function(test)) {
      keep <- test(data);      
    } else if (identical(modifier, "pattern")) {
      pattern <- as.character(test);
      keep <- (regexpr(pattern, data) != -1);  
    } else {
      keep <- (test == data);
    }
    df <- df[keep,,drop=FALSE];
  }

  if (!is.null(sortBy)) {
    args <- list();
    for (kk in seq(along=sortBy)) {
      name <- sortBy[kk];
      name <- strsplit(name, split=":")[[1]];
      if (length(name) > 1) {
        modifier <- name[1];
        name <- paste(name[-1], collapse=":");
      } else {
        modifier <- NULL;
      }
      cc <- pmatch(name, colnames(df));
      value <- df[,cc];
      if (!is.null(modifier)) {
        storage.mode(value) <- modifier;
      }
      args[[kk]] <- value;
    }
    o <- do.call("order", args=args);
    df <- df[o,];
  }

  df;
})

setMethodS3("getUnitIndices", "GdasAnnotationSet", function(this, ..., cdf=getCdf(this)) {
  df <- select(this, ...);
  names <- rownames(df);
  unitNames <- getUnitNames(cdf);
  pos <- match(names, unitNames);
  pos;
})

############################################################################
# HISTORY:
# 2006-08-28
# o Added select() [cool!] and then getUnitIndices().
# o Added getChipType() and getCdf().
# 2006-08-27
# o Created.
############################################################################
