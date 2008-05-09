setConstructorS3("GenericTabularFile", function(..., sep=c("\t", ","), quote="\"", fill=FALSE, skip=0, columnNames=TRUE, .verify=TRUE, verbose=FALSE) {
  # Argument 'columnNames':
  if (is.logical(columnNames)) {
    readColumnNames <- columnNames;
    columnNames <- NULL;
  } else if (is.character(columnNames)) {
    readColumnNames <- FALSE;
  } else {
    throw("Argument 'columnNames' must be either a logical or a character vector: ", class(columnNames)[1]);
  }

  this <- extend(GenericDataFile(...), "GenericTabularFile",
    .fileHeader = NULL,
    .columnNameTranslator = NULL,
    sep = sep,
    quote = quote,
    fill = fill,
    skip = skip,
    columnNames = columnNames,
    readColumnNames = readColumnNames
  );

  if (.verify)
    verify(this, ..., verbose=verbose);
  this;
})


setMethodS3("as.character", "GenericTabularFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", this, ...);
  class <- class(s);
  if (this$readColumnNames) {
    columns <- paste("'", getColumnNames(this), "'", sep="");
    s <- c(s, sprintf("Columns [%d]: %s", length(columns), paste(columns, collapse=", ")));
  } else {
    s <- c(s, sprintf("Columns [NA]: <not reading column names>"));
  }
  s <- c(s, sprintf("Number of data rows (lines): %d (%d)", nbrOfRows(this), nbrOfLines(this)));

  class(s) <- class;
  s;
})



setMethodS3("verify", "GenericTabularFile", function(this, ..., verbose=FALSE) {
  # Nothing to do?
  if (is.null(getPathname(this)))
    return(invisible(this));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Validating file contents");

  tryCatch({
    data <- readDataFrame(this, skip=this$skip, nrow=10, verbose=verbose);
  }, error = function(ex) {
    throw("File format error of the tabular file ('", getPathname(this), "'): ", ex$message);
  })

  verbose && exit(verbose);

  invisible(this);
}, private=TRUE)




setMethodS3("getColumnNameTranslator", "GenericTabularFile", function(this, ...) {
  this$.columnNameTranslator;
})


setMethodS3("setColumnNameTranslator", "GenericTabularFile", function(this, fcn, ...) {
  # Arguments 'fcn':
  if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", class(fcn)[1]);
  }

  this$.columnNameTranslator = fcn;
})


setMethodS3("translateColumnNames", "GenericTabularFile", function(this, names, ...) {
  nameTranslator <- getColumnNameTranslator(this);
  if (!is.null(nameTranslator)) {
    names <- nameTranslator(names);
    if (identical(attr(names, "isFinal"), TRUE))
      return(names);
  }

  # Do nothing
  names;
}, protected=TRUE)


setMethodS3("readColumnNames", "GenericTabularFile", function(this, ...) {
  as.logical(this$readColumnNames);
})


setMethodS3("getColumnNames", "GenericTabularFile", function(this, ...) {
  if (this$readColumnNames) {
    colnames <- getHeader(this, ...)$columns;
    colnames <- translateColumnNames(this, colnames);
  } else {
    colnames <- this$columnNames;
  }
  colnames;
})



setMethodS3("getHeader", "GenericTabularFile", function(this, ..., force=FALSE) {
  hdr <- this$.fileHeader;
  if (force || is.null(hdr)) {
    hdr <- readHeader(this, ...);
    this$.fileHeader <- hdr;
  }
  hdr;
})


setMethodS3("readHeader", "GenericTabularFile", function(this, con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading comment header & column header from ", class(this)[1]);

  # Open a file connection?
  if (is.null(con)) {
    pathname <- getPathname(this);
    verbose && cat(verbose, "Pathname: ", pathname);

    # Open file connection
    con <- file(pathname, open="r");
    on.exit({
      if (!is.null(con)) {
        close(con);
        con <- NULL;
      }
    })
  }


  ready <- FALSE;
  comments <- c();
  skip <- this$skip;
  while (!ready) {
    line <- readLines(con, n=1);
    isComments <- (regexpr("^#", line) != -1);
    if (!isComments) {
      if (skip == 0)
        break;
      skip <- skip - 1;
    }
    comments <- c(comments, line);
  }

  verbose && cat(verbose, "Header comments:", level=-20);
  verbose && str(verbose, comments, level=-20);

  # Infer column separator?
  sep <- this$sep;
  if (length(sep) > 1) {
    verbose && enter(verbose, "Identifying the separator that returns most columns");
    verbose && cat(verbose, "Separators:");
    verbose && str(verbose, sep);
    columns <- base::lapply(sep, FUN=function(split) {
      strsplit(line, split=split)[[1]];
    });
    nbrOfColumns <- sapply(columns, FUN=length);
    max <- which.max(nbrOfColumns);
    sep <- sep[max];
    verbose && printf(verbose, "Choosen separator: '%s' (0x%s)\n", sep, charToRaw(sep));
    verbose && exit(verbose);
  }

  if (this$readColumnNames) {
    verbose && print(verbose, line);
    columns <- strsplit(line, split=sep)[[1]];
    columns <- trim(columns);
    verbose && print(verbose, columns);
  } else {
    columns <- NULL;
  }

  # Remove quotes?
  quote <- this$quote;
  if (!is.null(quote)) {
    for (pattern in c(sprintf("^%s", quote), sprintf("%s$", quote))) {
      columns <- gsub(pattern, "", columns);
    }
  }

  verbose && cat(verbose, "Columns: ", paste(paste("'", columns, "'", sep=""), collapse=", "), level=-10);

  hdr <- list(
    comments=comments,
    columns=columns,
    sep=sep,
    quote=quote,
    skip=this$skip
  );

  verbose && str(verbose, hdr);

  verbose && exit(verbose);

  hdr;
}, protected=TRUE);



setMethodS3("getReadArguments", "GenericTabularFile", function(this, fileHeader=NULL, colClassPatterns=c("*"=NA), defColClass="NULL", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fileHeader':
  if (is.null(fileHeader)) {
    fileHeader <- getHeader(this);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Default column classes
  if (this$readColumnNames) {
    columns <- getColumnNames(this);
    nbrOfColumns <- length(columns);
    defColClasses <- rep(defColClass, nbrOfColumns);
    defColClassPatterns <- defColClasses;

    # Default columns?
    pos <- which(names(colClassPatterns) == "*");
    if (length(pos) > 0) {
      # Exclude extra '*':s
      if (length(pos) > 1) {
        colClassPatterns <- colClassPatterns[-(pos[-1])];
        pos <- pos[1];
      }
  
      # Insert defaults
      colClass <- colClassPatterns[pos];
      names <- names(colClassPatterns);
      if (length(colClassPatterns) > 1) {
        names <- insert(names[-pos], at=pos, values=rep("*", nbrOfColumns));
        idxs <- which(names == "*");
        names[idxs] <- sprintf("^%s$", columns);
  
        colClassPatterns <- insert(colClassPatterns[-pos], at=pos, 
                                   values=rep("*", nbrOfColumns));
        names(colClassPatterns) <- names;
        colClassPatterns[idxs] <- colClass;
      } else {
        colClassPatterns <- rep(colClass, nbrOfColumns);
        names(colClassPatterns) <- sprintf("^%s$", columns);
      }
    }

    verbose && cat(verbose, "Pattern used to identify column classes:", level=-20);
    verbose && print(verbose, colClassPatterns, level=-20);
  
    verbose && cat(verbose, "Generate column classes:");
    # Read everything by default
    colClasses <- defColClasses;
    names(colClasses) <- columns;
  
    # Update column classes according to patterns
    for (kk in seq(along=colClassPatterns)) {
      pattern <- names(colClassPatterns)[kk];
      idxs <- which(regexpr(pattern, columns) != -1);
      if (length(idxs) > 0) {
        colClass <- colClassPatterns[kk];
        colClasses[idxs] <- colClass;
      }
    }
  } else {
    colClasses <- NULL;
  }
  
  verbose && cat(verbose, "Column classes:", level=-20);
  verbose && print(verbose, colClasses, level=-20);

  # Inferred arguments
  args <- list(
    header      = FALSE,
    colClasses  = colClasses,
    sep         = fileHeader$sep,
    quote       = fileHeader$quote,
    fill        = this$fill,
    check.names = FALSE,
    na.strings  = c("---")
  );

  # Overwrite with user specified arguments, if any
  userArgs <- list(...);
  for (key in names(userArgs)) {
    args[[key]] <- userArgs[[key]];
  }

  args;
}, protected=TRUE);



setMethodS3("readDataFrame", "GenericTabularFile", function(this, con=NULL, rows=NULL, nrow=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rows':
  if (!is.null(rows)) {
    rows <- Arguments$getIndices(rows);
    nrow <- max(rows);
  }

  # Argument 'nrow':
  if (!is.null(nrow)) {
    nrow <- Arguments$getInteger(nrow, range=c(1,Inf));
  }
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading ", class(this)[1]);

  # Open a file connection?
  if (is.null(con)) {
    pathname <- getPathname(this);
    verbose && cat(verbose, "Pathname: ", pathname);

    # Open file connection
    con <- file(pathname, open="r");
    on.exit({
      if (!is.null(con)) {
        close(con);
        con <- NULL;
      }
    })
  }

  # Reading header
  hdr <- readHeader(this, con=con, verbose=verbose);

  # Get read arguments
  args <- getReadArguments(this, fileHeader=hdr, nrow=nrow, ..., verbose=less(verbose));
  args <- c(list(con), args);
  verbose && cat(verbose, "Arguments used to read tabular file:");
  verbose && str(verbose, args);

  columns <- getColumnNames(this);
  verbose && printf(verbose, "Column names (%d):\n", length(columns));
  verbose && cat(verbose, paste(columns, collapse=", "));

  if (!is.null(columns)) {
    verbose && enter(verbose, "Matching column names:");
    verbose && printf(verbose, "Column classes (%d):\n", length(args$colClasses));
    verbose && cat(verbose, paste(args$colClasses, collapse=", "));
    columns[args$colClasses == "NULL"] <- NA;
    columns <- columns[!is.na(columns)];
    verbose && exit(verbose);
  }

  # Read the data table
  verbose && enter(verbose, "Calling read.table()");
  verbose && print(verbose, args);
  data <- do.call("read.table", args=args);

  # Extract subset of rows?
  if (!is.null(rows)) {
    data <- data[rows,,drop=FALSE];
  } else {
    rownames(data) <- NULL;
  }

  if (length(columns) > 0) {
    if (ncol(data) != length(columns)) {
      throw("Number of read data columns does not match the number of column headers: ", ncol(data), " !=", length(columns));
    }
    colnames(data) <- columns;
  }

  verbose && str(verbose, data);
  verbose && exit(verbose);

  attr(data, "fileHeader") <- hdr;

  verbose && exit(verbose);

  data;
})

setMethodS3("readData", "GenericTabularFile", function(this, ...) {
  readDataFrame(this, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("nbrOfLines", "GenericTabularFile", function(this, ...) {
  pathname <- getPathname(this);
  con <- file(pathname, open="r");
  on.exit(close(con));
  count <- as.integer(0);
  while (TRUE) {
    dummy <- readLines(con=con, n=4096);
    n <- length(dummy);
    count <- count + n;
    if (n == 0)
      break;
  }
  count;
})


setMethodS3("readLines", "GenericTabularFile", function(con, ...) {
  # To please R CMD check
  this <- con;
  pathname <- getPathname(this);
  readLines(pathname, ...);
})


setMethodS3("nbrOfRows", "GenericTabularFile", function(this, ...) {
  hdr <- getHeader(this, ...);
  n <- nbrOfLines(this) - hdr$skip - as.integer(length(hdr$columns) > 0);
  n <- as.integer(n);
  n;
})


setMethodS3("nbrOfColumns", "GenericTabularFile", function(this, ...) {
  length(getColumnNames(this));
})


setMethodS3("dim", "GenericTabularFile", function(x) {
  # To please R CMD check.
  this <- x;

  c(nbrOfRows(this), nbrOfColumns(this));
})


############################################################################
# HISTORY:
# 2008-04-29
# o Added readLines(), nbrOfLines(), nbrOfRows() and dim().
# o Now readDataFrame() keeps the row names if arguments rows != NULL.
# 2008-04-25
# o Now argument 'verbose' of the constructor is passed to verfity().
# 2008-04-24
# o Added argument 'rows' to readDataFrame() for GenericTabularFile.
# 2008-04-14
# o Renamed readData() to readDataFrame() for GenericTabularFile.
# 2008-03-22
# o Added {get|set}ColumnNameTranslator().
# 2008-03-18
# o Now any '...' arguments to getReadArguments() override the inferred 
#   read arguments, e.g. na.strings="NA".
# 2008-02-27
# o Since 'affy' defines standardGeneric("colnames") and because S3 methods
#   are not found by such S4 generic functions, we avoid that method name,
#   and instead use getColumnNames().
# 2007-09-16
# o Removed all 'camelCaseNames' arguments.  Now column names are decided 
#   by getColumnNames() and translateColumnNames(), which can be overridden.
# 2007-09-14
# o Extracted from AffymetrixTabularFile.
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
