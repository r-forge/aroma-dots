setConstructorS3("IlluminaGenomeStudioTextFile", function(...) {
  extend(TabularTextFile(...), "IlluminaGenomeStudioTextFile");
})


setMethodS3("as.character", "IlluminaGenomeStudioTextFile", function(x, ...) {
  this <- x;
  s <- NextMethod(this, ...);
  class <- class(s);

  idx <- grep("^Columns", s);
  s <- s[seq(length=idx-1L)];

  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));

  nbrOfUnits <- NA;
  if (!is.null(this$.unitNames)) {
    nbrOfUnits <- nbrOfUnits(this);
  }
  s <- c(s, sprintf("Number of units: %d", nbrOfUnits));

  types <- getDataTypes(this);
  s <- c(s, sprintf("Data types [%d]: %s", length(types), hpaste(types)));

  sampleNames <- getSampleNames(this);
  s <- c(s, sprintf("Samples [%d]: %s", length(sampleNames), hpaste(sampleNames)));

  class(s) <- class;
  s;
})

setMethodS3("getPlatform", "IlluminaGenomeStudioTextFile", function(static, ...) {
  "Illumina";
}, static=TRUE)


setMethodS3("getDataSetName", "IlluminaGenomeStudioTextFile", function(this, ...) {
  path <- getPath(this);
  path <- getParent(path);
  basename(path);
})

setMethodS3("getChipType", "IlluminaGenomeStudioTextFile", function(this, ...) {
  path <- getPath(this);
  basename(path);
})


setMethodS3("getUnitNamesFile", "IlluminaGenomeStudioTextFile", function(this, force=FALSE, ...) {
  unf <- this$.unf;
  if (force || is.null(unf)) {
    chipType <- getChipType(this);
    nbrOfUnits <- nbrOfUnits(this);
    unf <- TextUnitNamesFile$byChipType(chipType, nbrOfUnits=nbrOfUnits, ...);
    this$.unf <- unf;
  }
  unf;
})


setMethodS3("getSampleNames", "IlluminaGenomeStudioTextFile", function(this, ...) {
  names <- getColumnNames(this, ...);
  pattern <- "[.](GType|Score|Log R Ratio|B Allele Freq)$";
  names <- grep(pattern, names, value=TRUE);
  names <- gsub(pattern, "", names);
  names <- unique(names);
  names;
})


setMethodS3("getDataTypes", "IlluminaGenomeStudioTextFile", function(this, toCamelCase=TRUE, ...) {
  sampleName <- getSampleNames(this)[1];
  pattern <- sprintf("^%s[.](.+)$", sampleName);
  names <- getColumnNames(this, ...);
  types <- grep(pattern, names, value=TRUE);
  types <- gsub(pattern, "\\1", types);
  keys <- toCamelCase(types);
  if (toCamelCase) {
    types <- keys;
  } else {
    names(types) <- keys;
  }
  types;
})


setMethodS3("nbrOfSamples", "IlluminaGenomeStudioTextFile", function(this, ...) {
  length(getSampleNames(this, ...));
})


setMethodS3("nbrOfUnits", "IlluminaGenomeStudioTextFile", function(this, ...) {
  length(getUnitNames(this, ...));
})

setMethodS3("getUnitNames", "IlluminaGenomeStudioTextFile", function(this, force=FALSE, ...) {
  unitNames <- this$.unitNames;

  if (force || is.null(unitNames)) {
    patterns <- c(Index="integer", Name="character");
    data <- readDataFrame(this, colClassPatterns=patterns);
    unitNames <- data$Name;
    this$.unitNames <- unitNames;
  }

  unitNames;
}) # getUnitNames()


setMethodS3("readUnitGenomePositions", "IlluminaGenomeStudioTextFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  patterns <- c(Index="integer", Name="character", Chr="character", Position="integer");
  data <- readDataFrame(this, colClassPatterns=patterns, ..., verbose=less(verbose, 10));
  unitNames <- data$Name;
  chrs <- data$Chr;
  chrs[chrs == "X"] <- 23L;
  chrs[chrs == "Y"] <- 24L;
  chrs[chrs == "M"] <- 25L;
  gp <- data.frame(chromosome=chrs, position=data$Position, stringsAsFactors=FALSE);
  rownames(gp) <- unitNames;
  gp;
}) # readUnitGenomePositions()



setMethodS3("readSampleData", "IlluminaGenomeStudioTextFile", function(this, sample, patterns=c("GType"="character", "Score"="double", "Log R Ratio"="double", "B Allele Freq"="double"), ..., verbose=FALSE) {
  # Argument 'sample':
  sampleNames <- getSampleNames(this);
  if (is.numeric(sample)) {
    sample <- Arguments$getIndex(sample, max=length(sampleNames));
    sample <- sampleNames[sample];
  } else if (!is.element(sample, sampleNames)) {
    throw("Unknown sample name: ", sample);    
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Read sample columns
  fields <- names(patterns);
  names(patterns) <- sprintf("^%s.%s", sample, names(patterns));
  patterns <- c(Name="character", patterns);

  data <- readDataFrame(this, colClassPatterns=patterns, ..., verbose=less(verbose, 5));
  unitNames <- data$Name;

  data <- data[-1L];
  rownames(data) <- unitNames;
  pattern <- sprintf("^%s.", sample);
  colnames(data) <- gsub(pattern, "", colnames(data));

  data;
})



setMethodS3("exportUnitNamesFile", "IlluminaGenomeStudioTextFile", function(this, chipType=getChipType(this), tags=NULL, filename=sprintf("%s,unitNames.txt", paste(c(chipType, tags), collapse=",")), path=file.path("annotationData", "chipTypes", chipType), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'tags':
  tags <- Arguments$getTags(tags);

  # Argument 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, ...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  unitNames <- getUnitNames(this);
  nbrOfUnits <- length(unitNames);

  hdr <- list(
    chipType=chipType,
    tags=tags,
    platform=getPlatform(this),
    nbrOfUnits=nbrOfUnits
  );

  data <- data.frame(unitName=unitNames);
  pathname <- writeDataFrame(data, file=pathname, path=NULL, ..., header=hdr);

  unf <- TextUnitNamesFile(pathname);

  invisible(unf);
})


setMethodS3("exportAromaUgpFile", "IlluminaGenomeStudioTextFile", function(this, unf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'unf':
  if (!is.null(unf)) {
    stopifnot(inherits(unf, "UnitNamesFile"));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Exporting genome position to UGP file");

  if (is.null(unf)) {
    unf <- getUnitNamesFile(this);
    verbose && cat(verbose, "Unit names file:");
    verbose && print(verbose, unf);
  }

  verbose && enter(verbose, "Reading genome positions");
  data <- readUnitGenomePositions(this, verbose=less(verbose, 10));
  verbose && str(verbose, data);
  verbose && exit(verbose);

  unitNames <- rownames(data);
  units <- indexOf(unf, names=unitNames);

  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);

  ugp <- AromaUgpFile$allocateFromUnitNamesFile(unf, ..., verbose=less(verbose, 10));
  verbose && print(verbose, ugp);

  ugp[units,1L] <- data$chromosome;
  ugp[units,2L] <- data$position;

  verbose && print(verbose, ugp);

  verbose && exit(verbose);

  invisible(ugp);
})


setMethodS3("extractTotalAndFracB", "IlluminaGenomeStudioTextFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading (total,fracB) data");

  patterns <- c("Log R Ratio"="double", "B Allele Freq"="double");
  data <- readSampleData(this, ..., patterns=patterns, verbose=less(verbose, 10));
  names <- colnames(data);
  names <- gsub("Log R Ratio", "total", names, fixed=TRUE);
  names <- gsub("B Allele Freq", "fracB", names, fixed=TRUE);
  colnames(data) <- names;

  # Transform
  data$total <- 2 * 2^data$total;

  data <- data[,c("total","fracB")];

  verbose && str(verbose, data);
  verbose && exit(verbose);
  data;
})


setMethodS3("exportTotalAndFracB", "IlluminaGenomeStudioTextFile", function(this, fields=c("total", "fracB"), samples=seq(length=nbrOfSamples(this)), unf=NULL, dataSet=getDataSetName(this), tags=NULL, rootPath="totalAndFracBData/", ..., overwrite=FALSE, force=FALSE, verbose=FALSE) {
  # Argument 'field':
  fields <- match.arg(fields, several.ok=TRUE);

  # Argument 'samples':
  samples <- Arguments$getIndices(samples, max=nbrOfSamples(this));

  # Argument 'unf':
  if (!is.null(unf)) {
    stopifnot(inherits(unf, "UnitNamesFile"));
  }

  # Arguments 'dataSet' & 'tags':
  dataSet <- Arguments$getCharacter(dataSet);
  tags <- Arguments$getTags(tags);

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

 
  units <- NULL;

  verbose && enter(verbose, "Exporting (total, fracB) data");
  verbose && cat(verbose, "Data set: ", dataSet);
  dataSetF <- paste(c(dataSet, tags), collapse=",");
  verbose && cat(verbose, "Fullname: ", dataSetF);
  verbose && cat(verbose, "Fields: ", paste(fields, collapse=", "));

  if (is.null(unf)) {
    unf <- getUnitNamesFile(this);
  }
  chipType <- getChipType(unf);
  chipType <- gsub(",unitNames", "", chipType, fixed=TRUE);
  verbose && cat(verbose, "Chip type: ", chipType);

  path <- file.path(rootPath, dataSetF, chipType);
  path <- Arguments$getWritablePath(path);
  verbose && cat(verbose, "Output path: ", path);

  sampleNames <- getSampleNames(this)[samples];
  nbrOfSamples <- length(sampleNames);
  verbose && cat(verbose, "Number of samples: ", nbrOfSamples);
  verbose && cat(verbose, "Samples:");
  verbose && str(verbose, sampleNames);

  verbose && cat(verbose, "Mapping to unit names:");
  verbose && print(verbose, unf);

  for (ii in seq(length=nbrOfSamples)) {
    sampleName <- sampleNames[ii];
    verbose && enter(verbose, sprintf("Sample %d ('%s') of %d", ii, sampleName, nbrOfSamples));

    fullname <- sampleName;

    footer <- list(
      srcFile=list(
        srcDataSet=dataSet,
        srcChipType=getChipType(this),
        srcFullName=getFullName(this),
        srcChecksum=getChecksum(this)
      )
    );

    data <- NULL;

    asbList <- list();
    for (field in fields) {
      # Identify output class
      if (field == "total") {
        signalClass <- AromaUnitTotalCnBinaryFile;
      } else if (field == "fracB") {
        signalClass <- AromaUnitFracBCnBinaryFile;
      }
  
      verbose && enter(verbose, "Exporting ", class(this)[1], " as an ", getName(signalClass));
      verbose && cat(verbose, "Signal: ", field);
  
      # Generate output filename
      filename <- sprintf("%s,%s.asb", fullname, field);
      pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE); 
      verbose && cat(verbose, "Output pathname: ", pathname);
      if (isFile(pathname)) {
        if (!overwrite) {
          verbose && cat(verbose, "Output file already exists. Return that instead.");
          asb <- signalClass$fromFile(pathname);
          asbList[[field]] <- asb;
          verbose && exit(verbose);
          next;
        }
      }
  
      verbose && cat(verbose, "File footer:");
      verbose && str(verbose, footer);

      # Reading data
      if (is.null(data)) {
        verbose && enter(verbose, "Reading data");
        data <- extractTotalAndFracB(this, sample=sampleName, verbose=less(verbose, 10));
        verbose && str(verbose, data);

        # Map to annotation data
        if (is.null(units)) {
          unitNames <- rownames(data);
          units <- indexOf(unf, names=unitNames);
          # Sanity check
          stopifnot(all(!is.na(units)));
        }
        verbose && exit(verbose);
      } # if (is.null(data)

      verbose && cat(verbose, "Values:");
      values <- data[,field, drop=TRUE];
      verbose && str(verbose, values);

      # Write to temporary file
      pathnameT <- pushTemporaryFile(pathname);

      verbose && enter(verbose, "Allocating output file"); 
      asb <- signalClass$allocateFromUnitNamesFile(unf, filename=pathnameT, path=NULL, footer=footer, ..., overwrite=overwrite, verbose=less(verbose, 25));
      verbose && print(verbose, asb);
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Writing data");
      asb[,1] <- values;
      verbose && exit(verbose);

      # Undo temporary filename
      pathname <- popTemporaryFile(pathnameT);

      asb <- signalClass(pathname);
      asbList[[field]] <- asb;

      rm(asb);
      verbose && exit(verbose);
    } # for (field ...)
    names(asbList) <- fields;
    rm(data);

    verbose && exit(verbose);
  } # for (ii ...)

  dsList <- list();
  for (field in fields) {
    # Identify output class
    if (field == "total") {
      signalClass <- AromaUnitTotalCnBinarySet;
      pattern <- ",total.asb$";
    } else if (field == "fracB") {
      signalClass <- AromaUnitFracBCnBinarySet;
      pattern <- ",fracB.asb$";
    }
    ds <- signalClass$byPath(path, pattern=pattern);
    idxs <- indexOf(ds, names=sampleNames);
    ds <- extract(ds, idxs);
    # Sanity check
    stopifnot(nbrOfFiles(ds) == length(sampleNames));
    dsList[[field]] <- ds;
  } # for (field ...)

  verbose && exit(verbose);

  invisible(dsList);
})


############################################################################
# HISTORY:
# 2012-02-17
# o Up for real usage.
# 2012-02-16
# o Now exporting to aroma binary files.
# 2012-02-13
# o Created.
############################################################################
