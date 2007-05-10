setConstructorS3("AptProbesetGenotype", function(dataSet=NULL, tags=c("APT", "BRLMM"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet")) {
      throw("Argument 'dataSet' is not an AffymetrixCelSet: ", 
                                                  class(dataSet)[1]);
    }

    nbrOfArrays <- nbrOfArrays(dataSet):
    if (nbrOfArrays < 6) {
      throw("Argument 'dataSet' contains too few samples. To fit the BRLMM model a minimum of six samples is required: ", nbrOfArrays);
    }
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }


  extend(Object(), "AptProbesetGenotype",
    .dataSet = dataSet,
    .tags = tags
  )
}, private=TRUE)


setMethodS3("as.character", "AptProbesetGenotype", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  s <- c(s, sprintf("Tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("convertTxtFilesToXdr", "AptProbesetGenotype", function(this, outPath=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'outPath':
  path <- getPath(this);
  path <- filePath(path, "txt");
  if (is.null(outPath)) {
    outPath <- getParent(getParent(path));
  } else {
    outPath <- Arguments$getWritablePath(outPath);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Converting TXT genotype call files into XDR files");
  verbose && cat(verbose, "Source path: ", path);
  verbose && cat(verbose, "Output path: ", outPath);

  # Assure data has been processed
  if (!isDone(this)) {
    throw("Data set not processed: ", getFullName(getDataSet(this)));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate all TXT genotype call files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating all TXT genotype call files");
  if (!isDirectory(path))
    throw("Path not found: ", path);
  pattern <- "[.]brlmm[.]txt$";
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);
  nbrOfFiles <- length(pathnames);
  verbose && cat(verbose, "Number of *.brlmm.txt files: ", nbrOfFiles);
  if (nbrOfFiles == 0)
    throw("Not *.brlmm.txt files found: ", path);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get CDF unit names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  unitNames <- getUnitNames(cdf);
  nbrOfUnits <- length(unitNames);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Convert each file to XDR
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outPathnames <- list(calls=c(), scores=c());
  for (kk in seq(along=pathnames)) {
    pathname <- pathnames[kk];

    # Get the sample name
    filename <- basename(pathname);
    pattern <- "[.]brlmm[.]txt$";
    fullname <- gsub(pattern, "", filename);
    name <- gsub(",.*$", "", fullname);

    verbose && enter(verbose, 
                 sprintf("Converting array %d ('%s') to XDR", kk, name));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Locate beginning of data table
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fh <- file(pathname, open="r");
    on.exit({
      if (!is.null(fh))
        close(fh);
    });
    skip <- 0;
    pattern <- "^SNPID";
    while(TRUE) {
      line <- readLines(con=fh, n=1);
      if (regexpr(pattern, line) != -1)
        break;
      skip <- skip + 1;
    }
    close(fh); fh <- NULL;


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Read data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading TXT file");
    verbose && cat(verbose, "Pathname: ", pathname);
    colClasses <- c("character", "integer", "double");
    data <- readTable(file=pathname, colClasses=colClasses, header=TRUE, 
                                                    sep="\t", skip=skip);
    verbose && exit(verbose);
    
    
    # Remap row to unit indices
    units <- match(data[,"SNPID"], unitNames);
    if (any(is.na(units))) {
      unknownUnits <- data[is.na(units),"SNPID"];
      throw("File contains unknown unit names: ", 
                                        paste(unitNames, collapse=", "));
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Genotype calls
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Calls");
    filename <- sprintf("%s,calls.xdr", name);
    pathname <- filePath(outPath, filename);
    verbose && cat(verbose, "Pathname: ", pathname);

    # Get calls
    values <- data$Call;  # NC=-1, AA=0, AB=1, BB=2

    # Convert to character strings
    labels <- c("-", "NC", "AA", "AB", "BB");
    values <- labels[values+3];

    # Expand to a complete CDF vector
    calls <- rep("-", nbrOfUnits);
    calls[units] <- values;
    rm(values);

    # Convert to factor
    calls <- as.factor(calls);

    attr(calls, "sampleName") <- name;  # Update the sample name

    verbose && str(verbose, calls);

    # Write to file
    saveObject(calls, file=pathname);

    outPathnames$calls <- c(outPathnames$calls, pathname);
    verbose && exit(verbose);



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Genotype scores
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Confidence scores");
    filename <- sprintf("%s,scores.xdr", name);
    pathname <- filePath(outPath, filename);
    verbose && cat(verbose, "Pathname: ", pathname);

    # Get complete CDF vector
    scores <- rep(NA, nbrOfUnits);
    scores[units] <- data$Score;

    attr(scores, "sampleName") <- name;  # Update the sample name

    verbose && str(verbose, scores);

    # Write to file
    saveObject(scores, file=pathname);

    outPathnames$scores <- c(outPathnames$scores, pathname);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...);

  verbose && exit(verbose);

  outPathnames;
}, private=TRUE)



setMethodS3("getGenotypeCallSet", "AptProbesetGenotype", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving genotype call set");

  path <- getPath(this);
  path <- getParent(path);
  callSet <- NULL;
  for (kk in 1:2) {
    tryCatch({
      callSet <- GenotypeCallXdrSet$fromFiles(path=path);
      return(callSet);
    }, error = function(ex) {
      # If XDR files are not available, try to convert from TXT files
      convertTxtFilesToXdr(this, outPath=path, verbose=less(verbose));
    })
  }

  if (is.null(callSet))
    throw("Failed to get GenotypeCallSet: ", path);

  verbose && exit(verbose);

  callSet;
})


setMethodS3("getRootPath", "AptProbesetGenotype", function(this, ...) {
  "genotypeData";
}, private=TRUE)


setMethodS3("getName", "AptProbesetGenotype", function(this, ...) {
  ds <- getDataSet(this);
  getName(ds);
})

setMethodS3("getFullName", "AptProbesetGenotype", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname;
})


setMethodS3("getPath", "AptProbesetGenotype", function(this, ...) {
  path <- filePath(getRootPath(this), getFullName(this), 
                getChipType(getCdf(this)), ".apt", expandLinks="any");
  path;
})

setMethodS3("getDataSet", "AptProbesetGenotype", function(this, ...) {
  this$.dataSet;
})

setMethodS3("nbrOfArrays", "AptProbesetGenotype", function(this, ...) {
  nbrOfArrays(getDataSet(this));
})

setMethodS3("getCdf", "AptProbesetGenotype", function(this, ...) {
  getCdf(getDataSet(this));
})

setMethodS3("getTags", "AptProbesetGenotype", function(this, collapse=NULL, ...) {
  ds <- getDataSet(this);
  tags <- c(getTags(ds), this$.tags);
  tags <- unique(tags);

  tags <- paste(tags, collapse=collapse);
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("getChrXFile", "AptProbesetGenotype", function(this, path=NULL, skip=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (!is.null(path))
    path <- Arguments$getWritablePath(path);


  # Get the CDF
  cdf <- getCdf(this);

  # Check if output file can be generated
  chipType <- getChipType(cdf);

  # Default output path
  if (is.null(path)) {
    path <- filePath("annotationData", "chipTypes", chipType, expandLinks="any");    
  }

  # Output pathname
  filename <- sprintf("%s.chrx", chipType);
  pathname <- filePath(path, filename);
  pathname <- Arguments$getWritablePathname(pathname);

  # Already done?
  if (skip && isFile(pathname)) {
    return(invisible(pathname));
  }

  # Retrieve the genome information
  gi <- getGenomeInformation(cdf);

  # Get the names of all units on chromosome X
  units <- getUnitsOnChromosome(gi, chromosome="X");
  units <- getUnitNames(cdf, units=units);

  # Exclude any NAs
  units <- units[!is.na(units)];

  # Write to file  
  lines <- c("all_chrx_no_par", units);
  writeLines(con=pathname, lines);

  # Return the pathname
  pathname;
}, private=TRUE)



setMethodS3("getPathnamesFile", "AptProbesetGenotype", function(this, ...) {
  # Get the data set
  ds <- getDataSet(this);

  # Create a temporary file
  path <- tempdir();
  # Not sure if it is ok to have a comma here. /HB 2006-12-20
  filename <- sprintf("%s_cel_file_list.txt", getName(ds));
  pathname <- file.path(path, filename);

  # Overwrite if already there
  if (isFile(pathname)) {
    file.remove(pathname);
    if (isFile(pathname))
      throw("Could not create temporary file: ", pathname);
  }

  # Write list of data file pathnames
  pathnames <- getPathnames(ds);
  lines <- c("cel_files", pathnames);
  writeLines(con=pathname, lines);

  # Return the pathname
  invisible(pathname);
}, private=TRUE)


setMethodS3("getAptPathnames", "AptProbesetGenotype", function(this, format=c("txt", "chp"), ...) {
  # Argument 'format':
  format <- match.arg(format);

  path <- getPath(this);
  path <- filePath(path, format);
  fullnames <- sapply(getDataSet(this), FUN=getFullName);
  filenames <- paste(fullnames, "brlmm", format, sep=".");
  pathnames <- file.path(path, filenames);
  pathnames;
}, private=TRUE);


setMethodS3("isDone", "AptProbesetGenotype", function(this, ...) {
  path <- getPath(this);
  path <- filePath(path, "chp");
  if (!isDirectory(path))
    return(FALSE);

  # Get full names of samples
  pathnames <- getAptPathnames(this, format="chp");
  isFile <- sapply(pathnames, FUN=isFile);

  all(isFile);
})


setMethodS3("readLog", "AptProbesetGenotype", function(this, collapse="\n", ...) {
  filename <- "apt-probeset-genotype.log";
  pathname <- filePath(getPath(this), filename);
  if (!isFile(pathname))
    throw("Log file does not exist: ", pathname);

  lines <- readLines(pathname);
  lines <- paste(lines, collapse=collapse);
  lines;
}, private=TRUE)


setMethodS3("showLog", "AptProbesetGenotype", function(this, ...) {
  log <- readLog(this, ...);
  displayCode(code=log);
})


setMethodS3("getExternalCommand", "AptProbesetGenotype", function(this, paths=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate Affymetrix Power Tools
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching for Affymetrix Power Tools");
  if (is.null(paths)) {
    paths <- c(".", 
               getOption("APT_PATH"), Sys.getenv("PATH_PATH"), 
               "apt/", "bin/");
    paths <- paste(paths, sep=";", collapse=";");
  }

  pattern <- "^apt-probeset-genotype";
  verbose && cat(verbose, "Filename pattern: ", pattern);

  pathname <- findFiles(pattern=pattern, paths=paths, firstOnly=TRUE);
  if (is.null(pathname)) 
    throw("Could not locate APT: ", pattern);
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && exit(verbose);

  # Try to call the APT command
  testCall <- sprintf("%s --version", pathname);
  res <- system(testCall);

  # Success?
  if (res != 0)
    throw("Failed to located/run APT binary: ", pathname);

  pathname;
}, private=TRUE)


setMethodS3("showAptVersion", "AptProbesetGenotype", function(this, ...) {
  res <- getExternalCommand(this, ...);
  invisible(res);
})



setMethodS3("process", "AptProbesetGenotype", function(this, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);


  verbose && enter(verbose, "Using APT to fit BRLMM");

  outPath <- getPath(this);
  verbose && cat(verbose, "Path: ", getPath(this));

  # Already done?
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Data already processed");
    callSet <- getGenotypeCallSet(this, verbose=verbose);
    return(callSet);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the APT command
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating external command");
  cmd <- getExternalCommand(this);
  verbose && cat(verbose, "Path: ", dirname(cmd));
  verbose && cat(verbose, "Command name: ", basename(cmd));
  outPath <- getPath(this);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create output path
  mkdirs(outPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Required options
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The pathname to the CDF file
  opts <- NULL;
  opts <- c(opts, sprintf("-c %s", getPathname(getCdf(this))));

  # The file listing all chromosome X units
  opts <- c(opts, sprintf("--chrX-snps %s", getChrXFile(this)));

  # The output path
  opts <- c(opts, sprintf("-o %s", outPath));

  # File listing the pathnames to all CEL files
  opts <- c(opts, sprintf("--cel-files %s", getPathnamesFile(this)));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Maximum verbosity
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opts <- c(opts, sprintf("--verbose 2"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The distribution used for quantile normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opts <- c(opts, sprintf("--write-sketch"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The chip-effect and probe-affinity esimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chip-effect estimates
  opts <- c(opts, sprintf("--summaries"));

  # Probe-affinity estimates
  opts <- c(opts, sprintf("--feat-effects"))

  # Residuals
#  opts <- c(opts, sprintf("--residuals"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Initial DM calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opts <- c(opts, sprintf("--dm-out"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Genotype calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opts <- c(opts, sprintf("--txt-output"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BRLMM parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Output the estimated priors
  opts <- c(opts, sprintf("--write-prior"));

  # Output the estimated cluster parameters for all SNPs
  opts <- c(opts, sprintf("--write-models"));


  callCmd <- paste(c(cmd, opts), collapse=" ");

  verbose && cat(verbose, "External call: ", callCmd);
  res <- system(callCmd);
  if (res != 0)
    throw("Error calling external command: ", callCmd);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return the GenotypeCallSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callSet <- getGenotypeCallSet(this, verbose=verbose);


  verbose && exit(verbose);

  callSet;
})


############################################################################
# HISTORY:
# 2007-01-06
# o Removed getChipType(). Use getCdf() first.
# 2006-12-20
# o For now, in case we will need it in the future, we output as much as
#   possible (except residuals) from APT.
# o Returned GenotypeCallSet is now compatible with what Crlmm is returning
#   in the sense that there is one *,calls.xdr and one *,scores.xdr file
#   per sample in both cases.
# o Default tags is not "APT,BRLMM".
# o Can now run APT and return a GenotypeCallSet.
# o Created.
############################################################################
