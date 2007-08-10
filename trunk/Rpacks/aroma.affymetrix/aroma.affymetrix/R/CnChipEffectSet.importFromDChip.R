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
#   \item{skip}{If @TRUE, already imported arrays will be skipped.}
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


  importUnits <- function(con, cdf, verbose=FALSE, ...) {
    verbose && enter(verbose, "Reading unit names (only)");
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
    verbose && enter(verbose, "Validating unit names towards the chip type");
    cdfUnitNames <- getUnitNames(cdf);
    units <- match(unitNames, cdfUnitNames);
    rm(cdfUnitNames);

    unknown <- which(is.na(units));
    nbrOfUnknown <- length(unknown);
    verbose && cat(verbose, "Number of unknown unit names: ", nbrOfUnknown);
    if (nbrOfUnknown == length(units)) {
      throw("None of the read unit names belongs to the '", getChipType(cdf),
                                                   "' CDF file: ", pathname);
    } 

    if (nbrOfUnknown > 0) {
      msg <- sprintf("Data file contains %d unknown unit names: %s", nbrOfUnknown, paste(unitNames[unknown], collapse=", "));
      throw(msg);
    }
    rm(unitNames, unknown);

    # Store only known units
    keep <- !is.na(units);
    units <- units[keep];
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    list(units=units, keep=keep);
  } # importUnits()


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
  monocellCdf <- clazz$createParamCdf(cdf);
  rm(cdf);
  verbose && print(verbose, monocellCdf);
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
    throw("dChip data file does not contain MBEI estimates: ", 
                                                        header$ModelMethod);
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
  nbrOfColumns <- length(colNames);
  verbose && printf(verbose, "Column names [%d]: %s\n", nbrOfColumns,
                                           paste(colNames, collapse=", "));

  # Infer sample names
  sampleNames <- gsub(" (call|SD)$", "", colNames[-1]);
  sampleNames <- unique(sampleNames);
  nbrOfSamples <- length(sampleNames);
  verbose && printf(verbose, "Sample names [%d]: %s\n", nbrOfSamples, paste(sampleNames, collapse=", "));


  nbrOfColsPerSample <- (nbrOfColumns-1)/nbrOfSamples;
  hasCalls <- any(regexpr(" call$", colNames[-1]) != -1);
  hasSDs <- any(regexpr(" SD$", colNames[-1]) != -1);
  verbose && printf(verbose, "Has calls: %s\n", hasCalls);
  verbose && printf(verbose, "Has SDs: %s\n", hasSDs);
  verbose && printf(verbose, "Columns per sample: %d\n", nbrOfColsPerSample);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Prepare to read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Record the current file position
  dataOffset <- seek(con, rw="read");
  verbose && cat(verbose, "Data file offset: ", dataOffset);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Import each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (signal, call, sd)
  sampleColClasses <- c("double", "NULL", "double")[1:nbrOfColsPerSample];
  cells <- NULL;

  verbose && enter(verbose, "Importing ", nbrOfSamples, " samples");
  for (kk in seq(length=nbrOfSamples)) {
    sampleName <- sampleNames[kk];
    verbose && enter(verbose, sprintf("Sample #%d (%s) of %d", 
                                            kk, sampleName, nbrOfSamples));

    # Create output filename
    filename <- sprintf("%s,chipEffects.CEL", sampleName);
    pathname <- file.path(outPath, filename);

    # Rename lower-case *.cel to *.CEL, if that is the case.  Old versions
    # of the package generated lower-case CEL files. /HB 2007-08-09
    pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

    verbose && cat(verbose, "Output pathname: ", pathname);

    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Already imported");
      verbose && exit(verbose);
      next;
    }

    cols <- 1 + nbrOfColsPerSample*(kk-1) + 1:nbrOfColsPerSample;
    verbose && cat(verbose, "Columns: ", seqToHumanReadable(cols));

    verbose && enter(verbose, "Retrieving chip-effect CEL file");
    verbose && cat(verbose, "Class: ", getName(clazz));
    if (isFile(pathname)) {
      cef <- clazz$fromFile(pathname, verbose=less(verbose));
    } else {
      cef <- clazz$fromDataFile(filename=filename, path=outPath, name=sampleName, cdf=monocellCdf, verbose=less(verbose));
    }
    cef$combineAlleles <- combineAlleles;
    cef$mergeStrands <- TRUE;
    verbose && print(verbose, cef);
    verbose && exit(verbose);

    if (is.null(cells)) {
      res <- importUnits(con, cdf=monocellCdf, verbose=verbose);
      units <- res$units;
      keep <- res$keep;
      rm(res);
      verbose && enter(verbose, "Getting CDF cell indices");
      cells <- getCellIndices(cef, units=units);
      verbose && cat(verbose, "CDF cell indices:");
      verbose && cat(verbose, "Number of units: ", length(cells));
      cells <- unlist(cells, use.names=FALSE);
      verbose && cat(verbose, "Number of cells: ", length(cells));
      verbose && str(verbose, cells);
      verbose && exit(verbose);
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && enter(verbose, "Reading data");
    colClasses <- rep("NULL", nbrOfColumns+1);
    colClasses[cols] <- sampleColClasses;
    seek(con, where=dataOffset, rw="read");
    data <- read.table(file=con, colClasses=colClasses, sep=sep, header=FALSE);
    data <- as.matrix(data[keep,,drop=FALSE]);
    dimnames(data) <- NULL;
    verbose && str(verbose, data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing chip effects");
    pathname <- getPathname(cef);
    rm(cef);

    if (hasSDs) {
      stdvs <- data[,2];
    } else {
      stdvs <- NULL;
    }
    updateCel(pathname, indices=cells, intensities=data[,1], stdvs=stdvs);
    rm(data, stdvs);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Define chip-effect set
  ces <- fromFiles(static, path=outPath, ...);
  setCdf(ces, monocellCdf);

  verbose && exit(verbose);

  ces;
}, static=TRUE, private=TRUE)


############################################################################
# HISTORY:
# 2007-08-09
# o CnChipEffectSet$importFromDChip() now creates CEL files with upper-case
#   filename extension "*.CEL", not "*.cel".  The reason for this is that
#   some software don't recognize lower case filename extensions :( 
# 2007-03-25
# o Made importFromDChip() more robust in identifying what columns goes with
#   what sample etc.
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
