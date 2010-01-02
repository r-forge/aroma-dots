###########################################################################/**
# @set "class=BroadBirdseedGenotypeTcgaDataFile"
# @RdocMethod exportGenotypeCallsAndConfidenceScores
#
# @title "Exports genotype calls and confidence scores"
#
#
# \description{
#    Export genotype calls and corresponding confidence scores for one sample
#    from a Birdseed data file to an 
#     @see "aroma.core::AromaUnitGenotypeCallFile" (genotype calls) and
#     an @see "aroma.core::AromaUnitSignalBinaryFile" (confidence scores).
#  }
#
# \arguments{
#   \item{dataSet}{@character, data set name.}
#   \item{unf}{A @see "aroma.core::UnitNamesFile".}
#   \item{...}{arguments passed to 
#      \code{AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(...)}.}
#   \item{rootPath}{@character, root path for output files.}
#   \item{force}{If @TRUE, existing output files are rewritten.}
#   \item{verbose}{@see "Verbose" object.}
# }
#
#
# \author{
#   Pierre Neuvial.
# }
#
#
# \details{
#   A @see "aroma.core::UnitNamesFile" is inferred from "chipType" using 
#   \code{AffymetrixCdfFile$byChipType()}.
#   Genotyping units (SNPs) are ordered according to this 
#   @see "aroma.core::UnitNamesFile".
#   Confidence scores are stored as one minus the input score as the
#   latter is a p-value.
# }
#
#*/###########################################################################
setMethodS3("exportGenotypeCallsAndConfidenceScores", "BroadBirdseedGenotypeTcgaDataFile", function(this, dataSet, unf, ..., rootPath="callData", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet, length=c(1,1));

  # Argument 'unf':
  unf <- Arguments$getInstanceOf(unf, "UnitNamesFile");

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

  fullname <- getFullName(this);
  verbose && cat(verbose, "Fullname: ", fullname);

  verbose && enter(verbose, "Exporting genotype calls and confidence scores for sample ", fullname);

  filename <- sprintf("%s,genotypes.acf", fullname);
  acfPathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE);
  if (!force && isFile(acfPathname)) {
    acf <- AromaUnitGenotypeCallFile(acfPathname);
    verbose && cat(verbose, "Already exported: ", getPathname(acf));
  } else {
    acf <- NULL;
  }

  filename <- sprintf("%s,confidenceScores.asb", fullname);
  asfPathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE);
  if (!force && isFile(asfPathname)) {
    asf <- AromaUnitSignalBinaryFile(asfPathname);
    verbose && cat(verbose, "Already exported: ", getPathname(asf));
  } else {
    asf <- NULL;
  }

  res <- list(acf=acf, asf=asf);

  # Nothing to do?
  if (!is.null(res$acf) && !is.null(res$asf)) {
    verbose && cat(verbose, "Nothing to do.");
    verbose && cat(verbose, "Existing data files:");
    verbose && print(verbose, res);
    verbose && exit(verbose);
    return(invisible(res));
  }



  verbose && enter(verbose, "Building file footer");
  srcFile <- list(filename = fullname,
                  fileSize = getFileSize(this),
                  checkSum = getChecksum(this));
  verbose && exit(verbose);
      
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transforming data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading data");
  verbose && enter(verbose, "Extracting calls and confidence scores");
  calls <- extractCalls(this);
  verbose && str(verbose, calls);
  verbose && print(verbose, summary.factor(calls));
    
  conf <- 1-extractConfidenceScores(this); # input data are p-values...
  verbose && str(verbose, conf);
  verbose && exit(verbose);
    
  verbose && enter(verbose, "Ordering according to UnitNamesFile");
  data <- readDataFrame(this);
  unitNames <- data[,"CompositeElement REF", drop=TRUE];
  rm(data);
  verbose && cat(verbose, "Unit names:");
  verbose && str(verbose, unitNames);
  units <- indexOf(unf, names=unitNames);
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  stopifnot(!anyMissing(units));
  verbose && exit(verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Step 1. Genotype calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(res$acf)) {
    verbose && enter(verbose, "Exporting genotype calls");
    pathname <- acfPathname;
    pathnameT <- sprintf("%s.tmp", pathname);
    pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);
    on.exit({
      if (isFile(pathnameT)) {
        file.remove(pathnameT);
      }
    }, add=TRUE);
      
    acf <- AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(unf, 
             filename=pathnameT, overwrite=force, ...);
      
    updateGenotypes(acf, units=units, calls=calls, encoding="birdseed", 
                    verbose=verbose);
    verbose && exit(verbose);
      
    footer <- readFooter(acf);
    footer$srcFile <- srcFile;
    writeFooter(acf, footer);
      
    # Rename temporary file
    file.rename(pathnameT, pathname);
    if (!isFile(pathname) || isFile(pathnameT)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }

    acf <- AromaUnitGenotypeCallFile(acfPathname);
    res$acf <- acf;
    verbose && exit(verbose);
  }
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Step 2. Confidence scores
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(res$asf)) {
    verbose && enter(verbose, "Exporting confidence scores");

    pathname <- asfPathname;
    pathnameT <- sprintf("%s.tmp", pathname);
    pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);
    on.exit({
      if (isFile(pathnameT)) {
        file.remove(pathnameT);
      }
    }, add=TRUE);
      
    asf <- AromaUnitSignalBinaryFile$allocateFromUnitNamesFile(unf, 
             filename=pathnameT, types="double", size=4, signed=TRUE, 
             overwrite=force, ...);
    naValue <- as.integer(NA);
    asf[, 1] <- naValue;
    asf[units, 1] <- conf;
    isNC <- (conf==0);
    if (sum(isNC)) {
      idxs <- units[whichVector(isNC)];
      asf[idxs, 1] <- naValue;
    }
    footer <- readFooter(asf);
    footer$srcFile <- srcFile;
    writeFooter(asf, footer);
      
    # Rename temporary file
    file.rename(pathnameT, pathname);
    if (!isFile(pathname) || isFile(pathnameT)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }

    asf <- AromaUnitSignalBinaryFile(asfPathname);
    res$asf <- asf;
    verbose && exit(verbose);
  } # if (...)

  verbose && cat(verbose, "Data files written:");
  verbose && print(verbose, res);
  
  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2009-11-02
# o CLEAN UP: Restructured code a bit.  Made less redundant.
# o Replaced argument 'chipType' with 'unf'; it's more generic.
# o BUG FIX: If the call file already exists, it was returned as a plain
#   AromaUnitSignalBinaryFile, not as an AromaUnitGenotypeCallFile.
# o BUG FIX: Called readDataFrame(df) instead of readDataFrame(this).
# o BUG FIX: Argument was named 'dataset' and not 'dataSet'.
# 2009-10-23
# o CLEANUPS.
# 2009-06-09
# o Genotype calls are now directly updated using:
#   'updateGenotypes(..., encoding="birdseed")' available in aroma.core v1.1.1.
# 2009-06-07
# o Added argument 'filenameFrom' to allow the output file name to be inferred 
#   either from the input file name or its header.
# o BUG FIX: NAs in input genotype files were converted to double
#   deletions (0,0). They are now converted to NCs ("No Call").
# 2009-05-28
# o A UnitNamesFile is now inferred from "chipType",
#   using "AffymetrixCdfFile$byChipType".
# 2009-05-19
# o Created.
############################################################################
