setConstructorS3("AromaAnnotationFile", function(...) {
  this <- extend(AffymetrixFile(...), "AromaAnnotationFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})

setMethodS3("as.character", "AromaAnnotationFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  s <- c(s, sprintf("Dimensions: %dx%d", nrow(this), ncol(this)));
  s <- c(s, sprintf("Column classes: %s", 
                         paste(getColClasses(this), collapse=", ")));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("readHeader", "AromaAnnotationFile", function(this, con=NULL, ...) {
  if (is.null(con)) {
    # Look for cached results
    hdr <- this$.hdr;
    if (!is.null(hdr))
      return(hdr);
  }

  knownDataTypes <- c("integer"=1, "double"=2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readBytes <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=1, n=n, signed=FALSE);
  }

  readShorts <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=2, n=n, signed=FALSE);
  }

  readInts <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=4, n=n, signed=FALSE);
  }

  readDataHeader <- function(con, ...) {
    # Number of elements (rows)
    nbrOfRows <- readInts(con);
    nbrOfRows <- Arguments$getInteger(nbrOfRows, range=c(0,100e6));
  
    # Number of fields (columns)
    nbrOfColumns <- readInts(con);
    nbrOfColumns <- Arguments$getInteger(nbrOfColumns, range=c(0,1000));

    # Types of columns
    types <- readBytes(con, n=nbrOfColumns);
    types <- Arguments$getIntegers(types, range=range(knownDataTypes));
    types <- names(knownDataTypes)[types];

    # Number of bytes per column
    sizes <- readBytes(con, n=nbrOfColumns);
    sizes <- Arguments$getIntegers(sizes, range=c(1,8));
    ok <- (sizes %in% c(1,2,4,8));
    if (any(!ok)) {
      cc <- which(!ok);
      throw("File format error. Detect one or more columns with invalid byte sizes, i.e. not in {1,2,4,8}: ", paste(paste(cc, sizes[cc], sep=":"), collapse=", "));
    }

    # Are the columns signed or not?
    signeds <- readBytes(con, n=nbrOfColumns);
    signeds <- Arguments$getIntegers(signeds, range=c(0,1));
    signeds <- as.logical(signeds);

    nbrOfBytes <- nbrOfRows*sizes;
    dataOffsets <- c(0, cumsum(nbrOfBytes[-length(nbrOfBytes)]));

    list(
      nbrOfRows=nbrOfRows, 
      nbrOfColumns=nbrOfColumns, 
      types=types, 
      sizes=sizes, 
      signeds=signeds,
      nbrOfBytes=nbrOfBytes,
      dataOffsets=dataOffsets,
      dataOffset=seek(con=con, rw="r")
    );
  } # readDataHeader()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open file?
  if (is.null(con)) {
    pathname <- getPathname(this);
    con <- file(pathname, open="rb");
    on.exit(close(con));
  }

  # Read magic
  trueMagic <- charToRaw("aroma");
  magic <- readBin(con=con, what=raw(), n=length(trueMagic));
  if (!identical(magic, trueMagic)) {
    asStr <- function(raw) {
      paste("[", paste(sprintf("%#0x", as.integer(magic)), collapse=","), "]", sep="");
    }
    throw("File format error. The read \"magic\" does not match the existing one: ", asStr(magic), " != ", asStr(trueMagic));
  }

  # File version
  version <- readBin(con=con, what=integer(), size=4, signed=FALSE);
  if (version < 0) {
    throw("File format error. Negative file version: ", version);
  }
  if (version > 10e3) {
    throw("File format error. Ridicolous large file version (>10e3): ", version);
  }

  if (version >= 1 && version <= 1) {
    dataHeader <- readDataHeader(con=con);
  } else {
    throw("Unknown file format version: ", version);
  }

  hdr <- list(dataHeader=dataHeader);

  # Cache result
  this$.hdr <- hdr;

  hdr;
}, protected=TRUE)



setMethodS3("readData", "AromaAnnotationFile", function(this, rows=NULL, columns=NULL, ..., verbose=FALSE) {
  # Open file
  pathname <- getPathname(this);
  con <- file(pathname, open="rb");
  on.exit(close(con));

  # Data header
  hdr <- readHeader(this, con=con)$dataHeader;
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'rows':
  if (is.null(rows)) {
    rows <- seq(length=hdr$nbrOfRows);
  } else if (is.logical(rows)) {
    rows <- which(rows);
    rows <- Arguments$getIndices(rows, range=c(1, hdr$nbrOfRows));
  } else {
    rows <- Arguments$getIndices(rows, range=c(1, hdr$nbrOfRows));
  }

  # Argument 'columns':
  if (is.null(columns)) {
    columns <- seq(length=hdr$nbrOfColumns);
  } else if (is.logical(columns)) {
    columns <- which(columns);
    columns <- Arguments$getIndices(columns, range=c(1, hdr$nbrOfColumns));
  } else {
    columns <- Arguments$getIndices(columns, range=c(1, hdr$nbrOfColumns));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Allocate return object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Allocating data object");
  colClasses <- hdr$types[columns];
  verbose && cat(verbose, "Number of rows: ", length(rows));
  verbose && cat(verbose, "Column classes: ", paste(colClasses, collapse=", "));
  data <- dataFrame(colClasses=colClasses, nrow=length(rows));
  colnames(data) <- make.unique(as.character(seq(length=ncol(data))));
  rownames(data) <- make.unique(as.character(rows));
  verbose && print(verbose, data, level=-30);
  verbose && exit(verbose);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading data");
  # First and last row to read in each column
  rrow <- range(rows);
  nbrOfRows <- as.integer(diff(rrow)+1);
  # Shift rows such that min(rows) == 1.
  rows <- rows - as.integer(rrow[1] - 1);

  # Record the current file offset
  dataOffsets <- hdr$dataOffsets[columns];

  # Read data in the order it appears on file
  o <- order(dataOffsets);

  count <- 0;
  for (kk in o) {
    count <- count + 1;
    verbose && enter(verbose, "Reading column #", count, " of ", length(o), level=-20);
    cc <- columns[kk];
    type <- hdr$types[cc];
    size <- hdr$sizes[cc];
    signed <- hdr$signeds[cc];

    verbose && printf(verbose, "Column %d: %s, %d bytes, signed=%s\n", cc, type, size, signed, level=-50);

    # Jump to the start of the data block
    dataOffset <- hdr$dataOffset + dataOffsets[kk] + (rrow[1]-1)*size;
    verbose && printf(verbose, "Data offset: %d\n", dataOffset, level=-50);
    seek(con=con, where=dataOffset, origin="start", rw="r");

    # Read from first to last row to be read, the discard unwanted.
    # TO DO: Optimize this.
    verbose && enter(verbose, "Reading binary data", level=-20);
    values <- readBin(con=con, n=nbrOfRows, what=type, size=size, 
                                         signed=signed, endian="little");
#    verbose && str(verbose, values, level=-30);
    verbose && exit(verbose);

    values <- values[rows];
#    verbose && str(verbose, values, level=-30);

    # Store data
    data[[o[kk]]] <- values;

    rm(values);
    gc <- gc();

    verbose && exit(verbose);
  }
  verbose && exit(verbose);
  
  data;
}, protected=TRUE)



setMethodS3("updateDataColumn", "AromaAnnotationFile", function(this, rows=NULL, column, values, .con=NULL, .hdr=NULL, .validateArgs=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  con <- .con;
  if (!is.null(con))
    seek(con, where=0, offset="start", rw="r"); 
  hdr <- .hdr;

  if (.validateArgs) {
    if (is.null(con)) {
      # Open file
      pathname <- getPathname(this);
      con <- file(pathname, open="r+b");
      on.exit(close(con));
    }
  
    # Data header
    if (is.null(hdr)) {
      hdr <- readHeader(this, con=con)$dataHeader;
    }

    # Argument 'rows':
    if (is.null(rows)) {
      rows <- seq(length=hdr$nbrOfRows);
    } else if (is.logical(rows)) {
      rows <- which(rows);
      rows <- Arguments$getIndices(rows, range=c(1, hdr$nbrOfRows));
    } else {
      rows <- Arguments$getIndices(rows, range=c(1, hdr$nbrOfRows));
    }
  
    # Argument 'column':
    column <- Arguments$getIndex(column, range=c(1, hdr$nbrOfColumns));

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }
  } # if (.validateArgs)

  verbose && enter(verbose, "Updating data column by writing to file");

  verbose && cat(verbose, "Number of rows: ", length(rows));
  verbose && cat(verbose, "Column: ", column);
  verbose && printf(verbose, "Values: %d %s(s)\n", length(values), mode(values));

  values <- rep(values, length.out=length(rows));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Prepare data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Optimizing data to be writing");
  verbose && cat(verbose, "Rows and values:");
  verbose && str(verbose, rows);
  verbose && str(verbose, values);

  # Remove duplicated rows
  rows <- rev(rows);
  values <- rev(values);
  dups <- duplicated(rows);
  rows <- rows[!dups];
  values <- values[!dups];
  rm(dups);

#  verbose && cat(verbose, "Rows and values:");
#  verbose && str(verbose, rows);
#  verbose && str(verbose, values);

  # Reorder rows
  o <- order(rows);
  rows <- rows[o];
  values <- values[o];
  rm(o);

  # Coerce data
  # Data type information
  type <- hdr$types[column];
  storage.mode(values) <- type;

  verbose && cat(verbose, "Rows and values:");
  verbose && str(verbose, rows);
  verbose && str(verbose, values);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Write data 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Shift rows such that min(rows) == 1.
  firstRow <- rows[1];
  rows <- rows - firstRow + 1;
  nbrOfRows <- rows[length(rows)];

  # Calculate the offset of the first element to read/write
  size <- hdr$sizes[column];
  dataOffset <- hdr$dataOffset + hdr$dataOffsets[column] + (firstRow-1)*size;

  # 1) Read existing data
  verbose && enter(verbose, "Reading existing data");
  seek(con=con, where=dataOffset, origin="start", rw="r");
  signed <- hdr$signeds[column];
  oldValues <- readBin(con=con, n=nbrOfRows, what=type, size=size, signed=signed, endian="little");
  verbose && str(verbose, oldValues);
  verbose && exit(verbose);
  
  # 2) Coerce and update the values
  storage.mode(oldValues) <- type;
  oldValues[rows] <- values;
  verbose && str(verbose, oldValues);
  rm(values, rows);

  # 3) Write back
  verbose && enter(verbose, "Writing updated data");
  seek(con=con, where=dataOffset, origin="start", rw="w");
  writeBin(con=con, object=oldValues, size=size, endian="little");
  flush(con);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(this);
})



setMethodS3("updateData", "AromaAnnotationFile", function(this, rows=NULL, columns=NULL, values, ..., verbose=FALSE) {
  # Open file
  pathname <- getPathname(this);
  con <- file(pathname, open="r+b");
  on.exit(close(con));

  # Data header
  hdr <- readHeader(this, con=con)$dataHeader;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'rows':
  if (is.null(rows)) {
    rows <- seq(length=hdr$nbrOfRows);
  } else if (is.logical(rows)) {
    rows <- which(rows);
    rows <- Arguments$getIndices(rows, range=c(1, hdr$nbrOfRows));
  } else {
    rows <- Arguments$getIndices(rows, range=c(1, hdr$nbrOfRows));
  }
  nbrOfRows <- length(rows);

  # Argument 'columns':
  if (is.null(columns)) {
    columns <- seq(length=hdr$nbrOfColumns);
  } else if (is.logical(columns)) {
    columns <- which(columns);
    columns <- Arguments$getIndices(columns, range=c(1, hdr$nbrOfColumns));
  } else {
    columns <- Arguments$getIndices(columns, range=c(1, hdr$nbrOfColumns));
  }
  nbrOfColumns <- length(columns);

  # Argument 'values':
  if (is.data.frame(values) || is.matrix(values)) {
    if (ncol(values) != nbrOfColumns) {
      throw("Number of columns in ", class(values), " 'values' does not match the number of specified columns: ", ncol(values), " != ", nbrOfColumns);
    }
  } else if (is.list(values)) {
    if (length(values) != nbrOfColumns) {
      throw("Number of elements in list 'values' does not match the number of specified columns: ", length(values), " != ", nbrOfColumns);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update each column
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the column in order, because that is faster
  o <- order(columns);
  count <- 0;
  for (kk in o) {
    count <- count + 1;
    verbose && enter(verbose, "Updating column #", count, " of ", length(o));
    cc <- o[kk];
    column <- columns[cc];

    # Extract the values
    if (is.data.frame(values) || is.matrix(values)) {
      theValues <- values[,cc];
    } else if (is.list(values)) {
      theValues <- values[[cc]];
    } else {
      # Is this strange?
      theValues <- values;
    }

    updateDataColumn(this, .con=con, .hdr=hdr, rows=rows, column=column, values=theValues, verbose=less(verbose));

    rm(theValues);

    verbose && exit(verbose);
  } # for (kk ...)

  invisible(this);
}, protected=TRUE)



setMethodS3("create", "AromaAnnotationFile", function(static, filename, path=NULL, nbrOfRows, types, sizes, signeds=FALSE, overwrite=FALSE, skip=FALSE, ...) {
  knownDataTypes <- c("integer"=1, "double"=2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  writeBytes <- function(con, values, ...) {
    values <- as.integer(values);
    writeBin(con=con, values, size=1, endian="little");
  }

  writeShorts <- function(con, values, ...) {
    values <- as.integer(values);
    writeBin(con=con, values, size=2, endian="little");
  }

  writeInts <- function(con, values, ...) {
    values <- as.integer(values);
    writeBin(con=con, values, size=4, endian="little");
  }

  writeFloats <- function(con, values, ...) {
    values <- as.double(values);
    writeBin(con=con, values, size=4, endian="little");
  }

  writeDoubles <- function(con, values, ...) {
    values <- as.double(values);
    writeBin(con=con, values, size=8, endian="little");
  }


  writeDataHeader <- function(con, nbrOfRows, types, sizes, signeds, ...) {
    # Number of elements (rows)
    writeInts(con=con, nbrOfRows);
  
    # Number of fields (columns)
    nbrOfColumns <- length(types);
    writeInts(con=con, nbrOfColumns);
  
    # Types of columns
    types <- knownDataTypes[types];
    writeBytes(con=con, types);

    # Number of bytes per column
    writeBytes(con=con, sizes);

    # Are the columns signed or not?
    writeBytes(con=con, signeds);
  } # writeDataHeader()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'nbrOfRows':
  nbrOfRows <- Arguments$getInteger(nbrOfRows, range=c(0,100e6));

  # Argument 'types':
  if (is.character(types))
    types <- knownDataTypes[types];
  types <- Arguments$getIntegers(types, range=range(knownDataTypes));
  types <- names(knownDataTypes)[types];

  nbrOfColumns <- length(types);
  nbrOfColumns <- Arguments$getInteger(nbrOfColumns, range=c(0,1000));

  # Argument 'sizes':
  sizes <- Arguments$getIntegers(sizes, range=c(1,8));
  ok <- (sizes %in% c(1,2,4,8));
  if (any(!ok)) {
    cc <- which(!ok);
    throw("Cannot create file. Detect one or more columns with invalid byte sizes, i.e. not in {1,2,4,8}: ", paste(paste(cc, sizes[cc], sep=":"), collapse=", "));
  }
  sizes <- rep(sizes, length.out=nbrOfColumns);

  # Argument 'signeds':
  signeds <- Arguments$getLogicals(signeds);
  signeds <- rep(signeds, length.out=nbrOfColumns);

  # Argument 'path':
  path <- Arguments$getWritablePath(path);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- Arguments$getWritablePathname(filename, path=path, 
                                         mustNotExist=(!overwrite && !skip));

  if (isFile(pathname)) {
    if (skip) {
      res <- newInstance(static, pathname);
      # TODO: We might retrieve an incompatible file.  Validate!
      return(res);
    } else if (!overwrite) {
      throw("Cannot not create file.  File already exists: ", pathname);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create empty file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open file
  con <- file(pathname, open="wb");
  on.exit(close(con));

  # Write magic
  magic <- charToRaw("aroma");
  writeBin(con=con, magic);

  # Write file version
  version <- 1;
  writeInts(con=con, version);

  # Write data header
  writeDataHeader(con=con, nbrOfRows=nbrOfRows, types=types, sizes=sizes, signeds=signeds);

  # Write empty data
  for (cc in seq(length=nbrOfColumns)) {
    size <- sizes[cc];
    type <- types[cc];
    value <- rep(NA, nbrOfRows);
    writeBin(con=con, value, size=size, endian="little");
    rm(value);
  }
  
  # Return
  newInstance(static, pathname);
}, static=TRUE)



setMethodS3("getColClasses", "AromaAnnotationFile", function(this, ...) {
  hdr <- readHeader(this)$dataHeader;
  hdr$types;
})


setMethodS3("dim", "AromaAnnotationFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  hdr <- readHeader(this)$dataHeader;
  c(hdr$nbrOfRows, hdr$nbrOfColumns);
})

setMethodS3("nrow", "AromaAnnotationFile", function(x, ...) {
  dim(x)[1];
})

setMethodS3("ncol", "AromaAnnotationFile", function(x, ...) {
  dim(x)[2];
})




setMethodS3("[", "AromaAnnotationFile", function(this, i=NULL, j=NULL, drop=FALSE) {
  # Read data
  data <- readData(this, rows=i, columns=j);

  # Drop dimensions?
  if (drop && length(dim(data)))
    data <- drop(data);
  
  data;
})


setMethodS3("[[", "AromaAnnotationFile", function(this, i) {
  if (!is.numeric(i))
    throw("Argument 'i' must be numeric: ", i);

  if (length(i) != 1)
    throw("Argument 'i' must be a single value: ", length(i));

  readData(this, columns=i)[[1]];
})



setMethodS3("[<-", "AromaAnnotationFile", function(this, i=NULL, j=NULL, value) {
  updateData(this, rows=i, columns=j, values=value, verbose=-20);
  invisible(this);
})


setMethodS3("subset", "AromaAnnotationFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  data <- readData(this);
  subset(data, ...);
})


setMethodS3("summary", "AromaAnnotationFile", function(object, ...) {
  # To please R CMD check
  this <- object;

  nbrOfColumns <- ncol(this);

  res <- lapply(seq(length=nbrOfColumns), FUN=function(cc) {
    s <- summary(this[,cc,drop=FALSE], ...);
    dimnames(s) <- NULL;
    s <- strsplit(s, split=":");
    names <- lapply(s, FUN=function(str) str[1]);
    values <- lapply(s, FUN=function(str) str[2]);
    names(values) <- names;
    values;
  })

  names <- sapply(res, FUN=function(s) names(s));
  unames <- unique(unlist(names, use.names=FALSE));
  emptyName <- paste(rep(" ", nchar(unames[1])+1), collapse="");

  for (kk in seq(along=res)) {
    s <- res[[kk]];
    emptyStr <- paste(rep(" ", nchar(s[[1]])), collapse="");
    thisNames <- names[[kk]];
    idx <- match(unames, thisNames);
    s <- s[idx];
    nok <- which(is.na(idx));
    s[nok] <- emptyStr;
    thisNames <- paste(thisNames, ":", sep="");
    thisNames[nok] <- emptyName;
    s <- paste(thisNames, s, sep="");
    res[[kk]] <- s;
  }

  res <- matrix(unlist(res, use.names=FALSE), ncol=nbrOfColumns);
  rownames(res) <- rep("", nrow(res));
  colnames(res) <- seq(length=ncol(res));
  class(res) <- "table";

  res;
})


setMethodS3("lapply", "AromaAnnotationFile", function(X, FUN, ...) {
  # To please R CMD check
  this <- X;

  nbrOfColumns <- ncol(this);
  res <- lapply(seq(length=nbrOfColumns), FUN=function(cc) {
    FUN(this[[cc]], ...);
  });

  res;
})


setMethodS3("colStats", "AromaAnnotationFile", function(this, FUN, ...) {
  res <- lapply(this, FUN=FUN, ...);
  res <- unlist(res, use.names=FALSE);
  res;
})

setMethodS3("colSums", "AromaAnnotationFile", function(x, ...) {
  colStats(x, FUN=sum, ...);
})

setMethodS3("colMeans", "AromaAnnotationFile", function(x, ...) {
  colStats(x, FUN=mean, ...);
})

setMethodS3("colMedians", "AromaAnnotationFile", function(x, ...) {
  colStats(x, FUN=median, ...);
})


############################################################################
# HISTORY:
# 2007-09-11
# o Added colStats(), colSums(), colMeans(), and colMedians().
# o Added dim(), nrow(), ncol() and a memory efficient lapply().
# o Created from AromaGenomeInformationFile.R.  This file type is more 
#   general.  There will be an option for a file header too, or rather a
#   file tail, because that is easier to expand.
############################################################################
