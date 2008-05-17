setConstructorS3("TabularTextFile", function(..., sep=c("\t", ","), quote="\"", fill=FALSE, skip=0, columnNames=TRUE, .verify=TRUE, verbose=FALSE) {
  # Argument 'columnNames':
  if (is.logical(columnNames)) {
    readColumnNames <- columnNames;
    columnNames <- NULL;
  } else if (is.character(columnNames)) {
    readColumnNames <- FALSE;
  } else {
    throw("Argument 'columnNames' must be either a logical or a character vector: ", class(columnNames)[1]);
  }

  this <- extend(GenericTabularFile(..., .verify=FALSE), "TabularTextFile",
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


setMethodS3("as.character", "TabularTextFile", function(x, ...) {
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



setMethodS3("verify", "TabularTextFile", function(this, ..., verbose=FALSE) {
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



setMethodS3("readColumnNames", "TabularTextFile", function(this, ...) {
  as.logical(this$readColumnNames);
})


setMethodS3("getColumnNames", "TabularTextFile", function(this, ...) {
  if (this$readColumnNames) {
    colnames <- getHeader(this, ...)$columns;
    colnames <- translateColumnNames(this, colnames);
  } else {
    colnames <- this$columnNames;
  }
  colnames;
})



setMethodS3("getHeader", "TabularTextFile", function(this, ..., force=FALSE) {
  hdr <- this$.fileHeader;
  if (force || is.null(hdr)) {
    hdr <- readHeader(this, ...);
    this$.fileHeader <- hdr;
  }
  hdr;
})


setMethodS3("readHeader", "TabularTextFile", function(this, con=NULL, ..., verbose=FALSE) {
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



setMethodS3("getReadArguments", "TabularTextFile", function(this, fileHeader=NULL, colClassPatterns=c("*"=NA), defColClass="NULL", ..., verbose=FALSE) {
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
  columns <- getColumnNames(this);
  if (!is.null(columns)) {
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
    args <- list(...);
    colClasses <- args$colClasses;
  }
  
  verbose && cat(verbose, "Column classes:", level=-20);
  verbose && print(verbose, colClasses, level=-20);

  # Inferred arguments
  args <- list(
    header      = this$readColumnNames,
    skip        = fileHeader$skip,
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



setMethodS3("readDataFrame", "TabularTextFile", function(this, con=NULL, rows=NULL, nrow=NULL, ..., verbose=FALSE) {
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading header to infer read.table() arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- readHeader(this, verbose=less(verbose, 5));

  # Get read arguments
  args <- getReadArguments(this, fileHeader=hdr, nrow=nrow, ..., 
                                               verbose=less(verbose, 5));

  verbose && cat(verbose, "Arguments inferred from file header:");
  verbose && print(verbose, args);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify names of columns read
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read the data using read.table()
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling read.table()");

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

  args <- c(list(con), args);
  verbose && cat(verbose, "Arguments used to read tabular file:");
  verbose && print(verbose, args);
  data <- do.call("read.table", args=args);
  verbose && cat(verbose, "Raw data read by read.table():");
  verbose && str(verbose, data);

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


setMethodS3("readColumns", "TabularTextFile", function(this, columns, colClasses=rep("character", length(columns)), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'columns':
  if (is.numeric(columns)) {
    columnNames <- getColumnNames(this);
    columns <- Arguments$getIndices(columns, range=c(1, length(columnNames)));
    columnNames <- columnNames[columns];
  } else {
    columnNames <- Arguments$getCharacters(columnNames);
  }


  colClassPatterns <- colClasses;
  names(colClassPatterns) <- sprintf("^%s$", columnNames);
  readDataFrame(this, colClassPatterns=colClassPatterns, ...);
}, protected=TRUE)



setMethodS3("nbrOfRows", "TabularTextFile", function(this, ...) {
  hdr <- getHeader(this, ...);
  n <- nbrOfLines(this) - hdr$skip - as.integer(length(hdr$columns) > 0);
  n <- as.integer(n);
  n;
})


setMethodS3("nbrOfLines", "TabularTextFile", function(this, ...) {
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


setMethodS3("readLines", "TabularTextFile", function(con, ...) {
  # To please R CMD check
  this <- con;
  pathname <- getPathname(this);
  readLines(pathname, ...);
})





############################################################################
# HISTORY:
# 2008-05-16
# o Took the text-file based parts of GenericTabularFile and placed them
#   in new subclass TabularTextFile.
# 2008-05-12
# o Added extractMatrix().
# o BUG FIX: getReadArguments() did not infer column classes if there was
#   no header to read but the column names was manually set.
# o BUG FIX: readDataFrame() did not read the first data row if there was
#   no column header; it was eaten up by a preceeding readHeader().
# 2008-04-29
# o Added readLines(), nbrOfLines(), nbrOfRows() and dim().
# o Now readDataFrame() keeps the row names if arguments rows != NULL.
# 2008-04-25
# o Now argument 'verbose' of the constructor is passed to verfity().
# 2008-04-24
# o Added argument 'rows' to readDataFrame() for TabularTextFile.
# 2008-04-14
# o Renamed readData() to readDataFrame() for TabularTextFile.
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
