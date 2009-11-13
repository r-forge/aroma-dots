setMethodS3("exportTotalAndFracB", "BroadTotalCopyNumberTcgaDataFile", function(this, dataSet, unf, ..., rootPath="totalAndFracBData", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'unf':
  if (!inherits(unf, "UnitNamesFile")) {
    throw("Argument 'unf' is not a UnitNamesFile: ", class(unf)[1]);
  }

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
  tags <- c("ratios", "total");

  # Unit indices (to be inferred)
  units <- NULL;

  # Export total signals
  allColumnNames <- getColumnNames(this);
  pattern <- ",Signal$";
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
    # Nothing to do?
    if (!force || isFile(pathname)) {
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
    data <- extractTotalCopyNumbers(this, drop=FALSE, 
                                            verbose=less(verbose,10));
    verbose && exit(verbose);

    # Map unit indices
    if (is.null(units)) {
      verbose && enter(verbose, "Mapping unit names to indices");
      unitNames <- rownames(data);
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
      column = cc,
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
# 2009-10-30
# o Added exportTotalAndFrac() for BroadTotalCopyNumberTcgaDataFile
#   and BroadTotalCopyNumberTcgaDataFileSet.
# o Created.
############################################################################
