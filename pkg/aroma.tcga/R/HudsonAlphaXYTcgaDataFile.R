setConstructorS3("HudsonAlphaXYTcgaDataFile", function(...) {
  this <- extend(TcgaDataFile(...), "HudsonAlphaXYTcgaDataFile");
  this;
})


setMethodS3("getExtensionPattern", "HudsonAlphaXYTcgaDataFile", function(static, ...) {
  "[.](XandYintensity[.]txt)$";
}, static=TRUE)



setMethodS3("getReadArguments", "HudsonAlphaXYTcgaDataFile", function(this, ..., colClassPatterns=c("*"="character", "(Chr|Pos)"="integer", "(X|Y)$"="double")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE)



setMethodS3("extractTheta", "HudsonAlphaXYTcgaDataFile", function(this, ..., drop=TRUE) {
  colClassPatterns <- c("CompositeElement REF"="character", "(X|Y)$"="double");

  throw("Not implemented");

  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);

  data;
})



setMethodS3("extractTotalAndFracB", "HudsonAlphaXYTcgaDataFile", function(this, ..., drop=TRUE) {
  data <- extractTheta(this, ..., drop=FALSE);

  # Add data[,2,] to the total data[,1,], if non-NA.
  ok <- whichVector(!is.na(data[,2,]));
  data[ok,1,] <- data[ok,1,] + data[ok,2,];
  rm(ok);

  # Calculate allele B fraction
  data[,2,] <- data[,2,] / data[,1,];

  # Drop singleton dimensions?
  nbrOfSamples <- dim[3];
  if (drop && nbrOfSamples == 1) {
    dimnames <- dimnames[-3];
    dim <- dim[-3];
  }

  dimnames(data)[[2]] <- c("total", "fracB");

  data;
})


setMethodS3("exportTotalAndFracB", "HudsonAlphaXYTcgaDataFile", function(this, dataSet, unf, ..., rootPath="totalAndFracBData", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet, length=c(1,1));

  # Argument 'unf':
  if (!inherits(unf, "UnitNamesFile")) {
    throw("Argument 'unf' is not a UnitNamesFile: ", class(unf)[1]);
  }

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Exporting ", class(this)[1]);
  chipType <- getChipType(unf, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);
  
  path <- file.path(rootPath, dataSet, chipType);
  path <- Arguments$getWritablePath(path);
  verbose && cat(verbose, "Path: ", path);

  # Unit indices (to be inferred)
  units <- NULL;

  allColumnNames <- getColumnNames(this);

  # Sample names to be imported
  pattern <- ",(X|Y)$";
  sampleNames <- grep(pattern, allColumnNames, value=TRUE);
  sampleNames <- gsub(pattern, "", sampleNames);
  # Sanity check
  stopifnot(length(sampleNames) %% 2 == 0);
  sampleNames <- unique(sampleNames);
  verbose && cat(verbose, "Samples to be processed:");
  verbose && print(verbose, sampleNames);
  nbrOfSamples <- length(sampleNames);

  typeList <- list(
    total = AromaUnitTotalCnBinarySet,
    fracB = AromaUnitFracBCnBinarySet
  );

  for (ii in seq(length=nbrOfSamples)) {
    sampleName <- sampleNames[ii];

    verbose && enter(verbose, sprintf("Exporting sample #%d ('%s') of %d", 
                                      ii, sampleName, nbrOfSamples));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Files to be exported
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    tags <- names(typeList);
    fullnames <- paste(sampleName, tags, sep=",");
    verbose && cat(verbose, "Full names to be exported:");
    verbose && print(verbose, fullnames);

    filenames <- sprintf("%s.asb", fullnames);
    pathnames <- file.path(path, filenames);

    res <- vector("list", length(typeList));
    names(res) <- names(typeList);
    for (kk in seq(along=typeList)) {
      name <- names(typeList)[kk];
      pathname <- pathnames[kk];
      verbose && cat(verbose, "Pathname: ", pathname);
      if (!force && isFile(pathname)) {
        verbose && cat(verbose, "Loading already existing file.");
        setClazz <- typeList[[kk]];
        className <- setClazz$getFileClass();
        fileClazz <- Class$forName(className);
        res[[name]] <- newInstance(fileClazz, pathname);
      }
    } # for (type ...)

    isDone <- (!sapply(res, FUN=is.null));
    if (!force && all(isDone)) {
      verbose && cat(verbose, "Sample already exported. Skipping.");
      verbose && exit(verbose);
      next;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Write to temporary files
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    pathnamesT <- sprintf("%s.tmp", pathnames);
    # Early sanity check
    for (pathnameT in pathnamesT) {
      if (isFile(pathnameT)) {
        throw("Temporary file already exists: ", pathnameT);
      }
    }
    

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Map unit indices
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if (is.null(units)) {
      verbose && enter(verbose, "Mapping unit names to  indices");
      colClassPatterns <- c("CompositeElement REF"="character");
      data <- readDataFrame(this, colClassPatterns=colClassPatterns);

      unitNames <- data[,1, drop=TRUE];
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


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Reading data for current sample
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    verbose && enter(verbose, "Reading data");
    verbose && cat(verbose, "Sample name: ", sampleName);
    pattern <- sprintf("%s,(X|Y)$", sampleName);
    verbose && cat(verbose, "Column pattern: ", pattern);
    colClassPatterns <- c("double");
    names(colClassPatterns) <- pattern;
    data <- readDataFrame(this, colClassPatterns=colClassPatterns);
    verbose && cat(verbose, "Data read:");
    verbose && str(verbose, data);
    columnNames <- colnames(data);
    verbose && cat(verbose, "Column names:");
    verbose && str(verbose, columnNames);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Generating tracability information
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    verbose && enter(verbose, "Generating 'srcFile' footer");
    srcFile <- list(
      filename = getFilename(this),
      filesize = getFileSize(this),
      checksum = getChecksum(this),
      sampleName = sampleName,
      columnNames = columnNames,
      valuesChecksum = digest(data)
    );
    verbose && str(verbose, srcFile);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Transforming data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    verbose && enter(verbose, "Transforming (X,Y) to (total,fracB)");
    data[,1] <- data[,1,drop=TRUE] + data[,2,drop=TRUE];
    data[,2] <- data[,2,drop=TRUE] / data[,1,drop=TRUE];
    colnames(data) <- names(typeList);
    verbose && str(verbose, data);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Writing data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    for (kk in seq(along=typeList)) {
      name <- names(typeList)[kk];

      verbose && enter(verbose, sprintf("Exporting #%d ('%s') of %d", 
                                          kk, name, length(typeList)));

      df <- res[[name]];
      if (!force && !is.null(df)) {
        verbose && cat(verbose, "Already exported. Skipping.");
        verbose && exit(verbose);
        next;
      }

      pathname <- pathnames[kk];
      verbose && cat(verbose, "Pathname: ", pathname);
      pathnameT <- pathnamesT[kk];
      verbose && cat(verbose, "Temporary pathname: ", pathnameT);
      # Sanity check
      if (!force && isFile(pathnameT)) {
        throw("Temporary file already exists: ", pathnameT);
      }

      on.exit({
        if (!is.null(pathnameT) && isFile(pathnameT)) {
          file.remove(pathnameT);
        }
      }, add=TRUE);


      setClazz <- typeList[[kk]];
      className <- setClazz$getFileClass();
      fileClazz <- Class$forName(className);

      verbose && enter(verbose, "Allocating temporary file");
      df <- fileClazz$allocateFromUnitNamesFile(unf, 
                                 filename=pathnameT, path=NULL);
      # Erase
      naValue <- as.double(NA);
      df[,1] <- naValue;
      footer <- readFooter(df);
      footer$srcFile <- srcFile;
      writeFooter(df, footer);
      verbose && exit(verbose);

      verbose && enter(verbose, "Writing signals");
      df[units,1] <- data[,name, drop=TRUE];
      verbose && exit(verbose);

      verbose && enter(verbose, "Renaming temporary file");
      # Rename temporary file
      file.rename(pathnameT, pathname);
      if (isFile(pathnameT) || !isFile(pathname)) {
        throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
      }
      pathnameT <- NULL;
      verbose && exit(verbose);

      df <- newInstance(fileClazz, pathname);
      verbose && cat(verbose, "Exported data file:");
      verbose && print(verbose, df);

      res[[name]] <- df;

      rm(df);
      verbose && exit(verbose);
    } # for (kk ...)
    rm(data);

    verbose && cat(verbose, "Sample name:", sampleName);
    verbose && cat(verbose, "Data files written:");
    verbose && print(verbose, res);

    rm(res);

    verbose && exit(verbose);
  } # for (ii ...)

  verbose && enter(verbose, "Setting up exported data sets");
  dsList <- list();
  for (kk in seq(along=typeList)) {
    name <- names(typeList)[kk];
    setClazz <- typeList[[kk]];
    pattern <- sprintf(",%s.asb$", tags[kk]);
    ds <- setClazz$byPath(path=path, pattern=pattern);
    dsList[[name]] <- ds;
  }
  verbose && print(verbose, dsList);
  verbose && exit(verbose);
  
  verbose && exit(verbose);

  invisible(dsList);
})


############################################################################
# HISTORY:
# 2009-11-06
# o Created.
############################################################################
