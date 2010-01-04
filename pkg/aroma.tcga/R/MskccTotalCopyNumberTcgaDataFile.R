setConstructorS3("MskccTotalCopyNumberTcgaDataFile", function(...) {
  this <- extend(TcgaDataFile(...), "MskccTotalCopyNumberTcgaDataFile");
  this;
})


setMethodS3("getReadArguments", "MskccTotalCopyNumberTcgaDataFile", function(this, ..., colClassPatterns=c("*"="character", "Pos$"="integer", "signal:Log2$"="double")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE)



setMethodS3("extractTotalLog2CopyNumbers", "MskccTotalCopyNumberTcgaDataFile", function(this, ..., drop=TRUE) {
  colClassPatterns <- c("ProbeID"="character", "signal:Log2$"="double");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);
  idx <- match("ProbeID", colnames(data));
  unitNames <- data[,idx];
  data <- data[,-idx,drop=FALSE];
  names <- names(data);

  pattern <- "(.*),(signal:Log2)$";
  sampleNames <- gsub(pattern, "\\1", names);
  nbrOfSamples <- length(sampleNames);

  # Coerce to a matrix  
  data <- as.matrix(data);
  rownames(data) <- unitNames;

  # A matrix? (probably never happens /HB 2009-08-23)
  if (drop && nbrOfSamples == 1) {
    data <- as.vector(data);
  }
  
  data;
})

setMethodS3("extractTotalCopyNumbers", "MskccTotalCopyNumberTcgaDataFile", function(this, ...) {
  data <- extractTotalLog2CopyNumbers(this, ...);

  # Transform to non-logged CNs
  data <- 2*2^data;

  if (is.matrix(data)) {
    names <- colnames(data);
    names <- gsub("signal:Log2", "ratio", names, fixed=TRUE);
    colnames(data) <- names;
  }

  data;
})


setMethodS3("exportTotal", "MskccTotalCopyNumberTcgaDataFile", function(this, dataSet, unf, ..., rootPath="totalAndFracBData", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'unf':
  unf <- Arguments$getInstanceOf(unf, "UnitNamesFile");

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Exporting ", class(this)[1]);

  chipType <- getChipType(unf, fullname=FALSE);
  
  path <- file.path(rootPath, dataSet, chipType);
  path <- Arguments$getWritablePath(path);

  # Tags added to each exported file data file
  tags <- c("log2ratio", "total");

  # Unit indices (to be inferred)
  units <- NULL;

  # Export total signals
  allColumnNames <- getColumnNames(this);
  pattern <- ",signal:Log2$";
  columnNames <- grep(pattern, allColumnNames, value=TRUE);
  verbose && cat(verbose, "Columns to be processed:");
  verbose && print(verbose, columnNames);

  for (cc in seq(along=columnNames)) {
    columnName <- columnNames[cc];

    verbose && enter(verbose, sprintf("Exporting column #%d ('%s') of %d", 
                     cc, columnName, length(columnNames)));

    # Sample name of exported data file
    sampleName <- gsub(pattern, "", columnName);

    fullname <- paste(c(sampleName, tags), collapse=",");
    verbose && cat(verbose, "Exported full name: ", fullname);

    filename <- sprintf("%s.asb", fullname);
    pathname <- file.path(path, filename);
    verbose && cat(verbose, "Export pathname: ", pathname);
    # Nothing to do?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Column already exported. Skipping.");
      verbose && exit(verbose);
      next;
    }
    # Export to a temporary file
    pathnameT <- sprintf("%s.tmp", pathname);
    if (isFile(pathnameT)) {
      throw("Temporary file already exists: ", pathnameT);
    }

    # Read data
    verbose && enter(verbose, "Reading column data");
    columnIdx <- match(columnName, allColumnNames);
    verbose && printf(verbose, "Column index: %d ('%s')\n", columnIdx, columnName);
    data <- readColumns(this, column=columnIdx, colClass="double")[,1];
    verbose && str(verbose, data);
    verbose && summary(verbose, data);
    verbose && exit(verbose);

    # Map unit indices
    if (is.null(units)) {
      verbose && enter(verbose, "Mapping unit names to indices");
      colClassPatterns <- c("ProbeID"="character");
      unitNames <- readDataFrame(this, colClassPatterns=colClassPatterns)[,1];
      verbose && cat(verbose, "Unit names:");
      verbose && str(verbose, unitNames);

      verbose && print(verbose, unf);
      units <- indexOf(unf, names=unitNames);
      verbose && cat(verbose, "Units:");
      verbose && str(verbose, units);

      # Sanity check
      if (anyMissing(units)) {
        missing <- unitNames[is.na(units)];
        throw("There are ", length(missing), " unknown unit names: ", 
                                 paste(head(missing, 3), collapse=", "));
      }
      verbose && exit(verbose);
    }

    # Drop attributes
    data <- as.vector(data);

    verbose && enter(verbose, "Generating 'srcFile' footer");
    srcFile <- list(
      filename = getFilename(this),
      filesize = getFileSize(this),
      checksum = getChecksum(this),
      column = columnIdx,
      columnName = columnName,
      valuesChecksum = digest(data)
    );
    verbose && str(verbose, srcFile);
    verbose && exit(verbose);

    on.exit({
      if (!is.null(pathnameT) && isFile(pathnameT)) {
        file.remove(pathnameT);
      }
    }, add=TRUE);

    verbose && enter(verbose, "Allocating temporary file");
    df <- AromaUnitTotalCnBinaryFile$allocateFromUnitNamesFile(unf, 
                                            filename=pathnameT, path=NULL);
    footer <- readFooter(df);
    footer$srcFile <- srcFile;
    writeFooter(df, footer);
    verbose && exit(verbose);

    verbose && enter(verbose, "Write signals");
    df[units,1] <- data;
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Renaming temporary file");
    # Rename temporary file
    file.rename(pathnameT, pathname);
    if (isFile(pathnameT) || !isFile(pathname)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    pathnameT <- NULL;
    verbose && exit(verbose);

    # Validate
    df <- AromaUnitTotalCnBinaryFile(pathname);
    verbose && cat(verbose, "Exported data file:");
    verbose && print(verbose, df);

    verbose && exit(verbose);
  } # for (cc ...)

  ds <- AromaUnitTotalCnBinarySet$byPath(path);
  verbose && print(verbose, ds);

  verbose && exit(verbose);

  invisible(ds);
})


############################################################################
# HISTORY:
# 2010-01-03
# o Created from HarvardTotalCopyNumberTcgaDataFile.R.
############################################################################
