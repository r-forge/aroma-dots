setConstructorS3("AffymetrixTabularFile", function(..., sep=c("\t", ","), quote="\"", fill=FALSE, .verify=TRUE) {
  this <- extend(AffymetrixFile(...), "AffymetrixTabularFile",
    .header = NULL,
    sep = sep,
    quote = quote,
    fill = fill
  );

  if (.verify)
    verify(this, ...);
  this;
})


setMethodS3("as.character", "AffymetrixTabularFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", this, ...);
  class <- class(s);
  header <- getHeader(this);
  columns <- paste("'", header$columns, "'", sep="");
  s <- c(s, sprintf("Columns [%d]: %s", length(columns), paste(columns, collapse=", ")));
  class(s) <- class;
  s;
})


setMethodS3("findByChipType", "AffymetrixTabularFile", function(static, chipType, pattern=sprintf("^%s.*[.]...$", chipType), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findAnnotationDataByChipType(chipType, pattern);
  pathname;
}, static=TRUE, protected=TRUE)



setMethodS3("fromChipType", "AffymetrixTabularFile", function(static, chipType, ...) {
  # Search for the genome information file
  pathname <- static$findByChipType(chipType, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix tabular file: ", chipType);
  newInstance(static, pathname);
}, static=TRUE)


setMethodS3("verify", "AffymetrixTabularFile", function(this, ..., verbose=FALSE) {
  # Nothing to do?
  if (is.null(getPathname(this)))
    return(invisible(this));

  tryCatch({
    data <- readData(this, nrow=10, verbose=verbose);
  }, error = function(ex) {
    throw("File format error of the Affymetrix tabular file: ", getPathname(this));
  })

  invisible(this);
}, private=TRUE)


setMethodS3("colnames", "AffymetrixTabularFile", function(this, ...) {
  getHeader(this, ...)$columns;
})


setMethodS3("getHeader", "AffymetrixTabularFile", function(this, ..., force=FALSE) {
  header <- this$.header;
  if (force || is.null(header)) {
    header <- readHeader(this, ...);
    this$.header <- header;
  }
  header;
})

setMethodS3("readHeader", "AffymetrixTabularFile", function(this, con=NULL, ..., verbose=FALSE) {
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
  while (!ready) {
    line <- readLines(con, n=1);
    isComments <- (regexpr("^#", line) != -1);
    if (!isComments)
      break;
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
    columns <- lapply(sep, FUN=function(split) {
      strsplit(line, split=split)[[1]];
    });
    nbrOfColumns <- sapply(columns, FUN=length);
    max <- which.max(nbrOfColumns);
    sep <- sep[max];
    verbose && printf(verbose, "Choosen separator: '%s' (0x%s)\n", sep, charToRaw(sep));
    verbose && exit(verbose);
  }

  columns <- strsplit(line, split=sep)[[1]];
  columns <- trim(columns);

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
    quote=quote
  );

  verbose && str(verbose, hdr);

  verbose && exit(verbose);

  hdr;
}, protected=TRUE);


setMethodS3("getReadArguments", "AffymetrixTabularFile", function(this, header, colClassPatterns=c("*"=NA), defColClass="NULL", camelCaseNames=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  columns <- header$columns;
  if (camelCaseNames)
    columns <- toCamelCase(columns);

  # Default column classes
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

  verbose && cat(verbose, "Column classes:", level=-20);
  verbose && print(verbose, colClasses, level=-20);

  args <- list(
    header      = FALSE,
    colClasses  = colClasses,
    sep         = header$sep,
    quote       = header$quote,
    fill        = this$fill,
    check.names = FALSE,
    na.strings  = c("---"),
    ...
  );

  args;
}, protected=TRUE);


setMethodS3("readData", "AffymetrixTabularFile", function(this, con=NULL, camelCaseNames=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  args <- getReadArguments(this, header=hdr, camelCaseNames=camelCaseNames, ..., verbose=less(verbose));
  args <- c(list(con), args);
  verbose && cat(verbose, "Arguments used to read tabular file:");
  verbose && str(verbose, args);

  columns <- hdr$columns;
  columns[args$colClasses == "NULL"] <- NA;
  columns <- columns[!is.na(columns)];

  # Make column names into camelCase names?
  if (camelCaseNames) {
    columns <- toCamelCase(columns);
  }

  # Read the data table
  verbose && enter(verbose, "Calling read.table()");
  data <- do.call("read.table", args=args);
  rownames(data) <- NULL;
  if (ncol(data) != length(columns)) {
    throw("Number of read data columns does not match the number of column headers: ", ncol(data), " !=", length(columns));
  }
  colnames(data) <- columns;
  verbose && str(verbose, data);
  verbose && exit(verbose);

  attr(data, "header") <- hdr;

  verbose && exit(verbose);

  data;
})


############################################################################
# HISTORY:
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
