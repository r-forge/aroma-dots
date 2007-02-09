###########################################################################/**
# @set "class=CnChipEffectSet"
# @RdocMethod importFromDChip
#
# @title "Imports copy-number estimates from a dChip result file"
#
# \description{
#  @get "title".
#  Currently only total copy-number estimates can be imported, that is
#  if dChip fitted the PLM with allele A and allele B combined.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the dChip result file.}
#   \item{path}{An optional path to the file.}
#   \item{combineAlleles}{If @TRUE, .}
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
#   \item{skip}{If @TRUE, already imported chip effects will not be imported
#     again.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "CnChipEffectSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("importFromDChip", "CnChipEffectSet", function(static, filename, path=NULL, combineAlleles=TRUE, cdf, ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readMbeiHeader <- function(con, ...) {
    res <- list();
    while(TRUE) {
      line <- readLines(con, n=1);
      if (line == "[Data]")
        break;
      line <- strsplit(line, split="=")[[1]];
      name <- line[1];
      value <- paste(line[-1], collapse="=");
      res[[name]] <- value;
    }
    
    res;    
  } # readMbeiHeader()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Argument 'combineAlleles':
  combineAlleles <- Arguments$getLogical(combineAlleles);
  if (!combineAlleles)
    throw("Currently only 'combineAlleles=TRUE' is supported");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  ces <- NULL;

  verbose && enter(verbose, "Importing MBEI estimates from dChip tabular file");
  verbose && cat(verbose, "Pathname: ", pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Infer data path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rootPath <- "plmData";
  chipType <- getChipType(cdf);
  currPath <- dirname(pathname);
  while(TRUE) {
    dataSetName <- basename(currPath);
    if (regexpr(chipType, dataSetName) == -1)
      break;
    currPath <- dirname(currPath);
  }

  outPath <- file.path(rootPath, dataSetName, chipType);
  outPath <- Arguments$getWritablePath(outPath);

  verbose && cat(verbose, "Output path: ", outPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get CDF for chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving CDF for chip effects");
  verbose && printf(verbose, "Chip type: %s,monocell\n", chipType);
   # Get the ChipEffectFile class specific for this set
  clazz <- getChipEffectFileClass(static);
  ceCdf <- clazz$createParamCdf(cdf);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open file & assert file format
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  con <- file(pathname, open="r");
  on.exit(close(con));

  line <- readLines(con, n=1);
  magic <- "[dChip Expression Data File]";
  if (line != magic) {
    throw(sprintf("File format error: First line is not '%s': %s", 
                                                             magic, line));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read and validate header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  header <- readMbeiHeader(con);
  verbose && cat(verbose, "File header:");
  verbose && str(verbose, header);

  if (tolower(header$ModelMethod) != tolower("Model-based expression")) {
    throw("dChip data file does not contain MBEI estimates: ", header$ModelMethod);
  }

  isLog <- (regexpr("^no", tolower(header$LogTransformed)) == -1);
  verbose && cat(verbose, "Log-transformed: ", isLog);
  if (isLog) {
    throw("dChip data file contains log-transformed values. Not supported.");
  }

  hasStddev <- (tolower(header$OutputBothSignalAndCall) == "yes");
  verbose && cat(verbose, "Has standard-deviation estimates: ", hasStddev);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read column names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  colNames <- readLines(con, n=1);

  # Guess column separator
  if (regexpr("\t", colNames) != -1) {
    sep <- "\t";
    verbose && cat(verbose, "File is tab-delimited");
  } else {
    sep <- ",";
    verbose && cat(verbose, "File is comma-separated");
  }

  colNames <- unlist(strsplit(colNames, split=sep));
  verbose && cat(verbose, "Column names: ", paste(colNames, collapse=", "));

  nbrOfColsPerSample <- ifelse(hasStddev, 3, 2);
  sampleNames <- colNames[seq(from=2, to=length(colNames), by=nbrOfColsPerSample)];
  nbrOfSamples <- length(sampleNames);
  verbose && cat(verbose, "Number of samples: ", nbrOfSamples);
  verbose && cat(verbose, "Sample names: ", paste(sampleNames, collapse=", "));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Prepare to read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfColumns <- length(colNames);
  colClasses <- rep("double", nbrOfColumns);
  # First column contains probeset names
  colClasses[1] <- "character";
  # Call column contains strings
  colClasses[grep("call", colNames)] <- "character";
  verbose && cat(verbose, "Column classes: ", paste(colClasses, collapse=", "));
  # Skip the last empty column (due to the extra tab outputted by dChip)
  colClasses <- c(colClasses, "NULL");

  # Record the current file position
  dataOffset <- seek(con, rw="read");
  verbose && cat(verbose, "Data file offset: ", dataOffset);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read unit names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading unit names");
  seek(con, where=dataOffset, rw="read");
  # Skip the last empty column (due to the extra tab outputted by dChip)
  colClasses <- rep("NULL", nbrOfColumns+1);
  colClasses[1] <- "character";
  unitNames <- read.table(file=con, colClasses=colClasses, sep=sep, header=FALSE);
  unitNames <- unlist(unitNames, use.names=FALSE);

  nbrOfUnits <- length(unitNames);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  verbose && cat(verbose, "Unit names: ");
  verbose && str(verbose, unitNames);
  verbose && exit(verbose);

  # Get the dChip-to-CDF unit map
  cdfUnitNames <- getUnitNames(cdf);
  units <- match(unitNames, cdfUnitNames);
  unknown <- which(is.na(units));
  nbrOfUnknown <- length(unknown);
  verbose && cat(verbose, "Number of unknown unit names: ", nbrOfUnknown);
  if (nbrOfUnknown == length(units)) {
    throw("Non of the read unit names belongs to the '", chipType, "' CDF file: ", pathname);
  } 

  if (nbrOfUnknown > 0) {
    msg <- sprintf("Data file contains %d unknown unit names: %s", nbrOfUnknown, paste(unitNames[unknown], collapse=", "));
    throw(msg);
  }

  # Store only known units
  keep <- !is.na(units);
  units <- units[keep];
  nbrOfUnits <- length(units);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Import each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sampleColClasses <- c("double", "NULL", "double")[1:nbrOfColsPerSample];
  cells <- NULL;

  verbose && enter(verbose, "Importing ", nbrOfSamples, " samples");
  for (kk in seq(length=nbrOfSamples)) {
    sampleName <- sampleNames[kk];
    verbose && enter(verbose, sprintf("Sample #%d (%s)", kk, sampleName));

    # Create output filename
    filename <- sprintf("%s,chipEffects.cel", sampleName);
    pathname <- file.path(outPath, filename);
    verbose && cat(verbose, "Output pathname: ", pathname);

    if (isFile(pathname) && skip) {
      verbose && cat(verbose, "Already imported");
      verbose && exit(verbose);
      next;
    }

    cols <- 1 + nbrOfColsPerSample*(kk-1) + 1:nbrOfColsPerSample;
    verbose && cat(verbose, "Columns: ", seqToHumanReadable(cols));

    verbose && enter(verbose, "Retrieving chip-effect CEL file");
    if (isFile(pathname)) {
      cef <- clazz$fromFile(pathname, verbose=less(verbose));
    } else {
      cef <- clazz$fromDataFile(filename=filename, path=outPath, name=sampleName, cdf=ceCdf, verbose=less(verbose));
    }
    cef$combineAlleles <- combineAlleles;
    cef$mergeStrands <- TRUE;
    verbose && print(verbose, cef);
    verbose && exit(verbose);

    if (is.null(cells)) {
      verbose && enter(verbose, "Getting CDF cell indices");
      cells <- getCellIndices(cef, units=units);
      verbose && cat(verbose, "CDF cell indices:");
      verbose && cat(verbose, "Number of units: ", length(cells));
      cells <- unlist(cells, use.names=FALSE);
      verbose && cat(verbose, "Number of cells: ", length(cells));
      verbose && str(verbose, cells);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Reading");
    colClasses <- rep("NULL", nbrOfColumns+1);
    colClasses[cols] <- sampleColClasses;
    seek(con, where=dataOffset, rw="read");
    data <- read.table(file=con, colClasses=colClasses, sep=sep, header=FALSE);
    data <- as.matrix(data[keep,]);
    dimnames(data) <- NULL;
    verbose && str(verbose, data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing chip effects");
    pathname <- getPathname(cef);
    updateCel(pathname, indices=cells, intensities=data[,1], stdvs=data[,2]);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Define chip-effect set
  ces <- fromFiles(static, path=outPath);
  ces$combineAlleles <- combineAlleles;
  ces$mergeStrands <- TRUE;

  verbose && exit(verbose);

  ces;
}, static=TRUE, private=TRUE)


############################################################################
# HISTORY:
# 2007-02-03
# o Added Rdoc comments.
# 2007-01-03
# o Now, if 'skip=FALSE' and chip-effect file already exists, a new file
#   is not created, but instead its contents is updated.
# 2007-01-02
# o Verified to work with SNP data modelled as PM=PMA+PMB (combined alleles)
#   exported from dChip v2006-12-14.  Imports both chip effects and standard
#   deviations of such, i.e. 'theta' and 'sdTheta'.
# o Created.
############################################################################
