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
#   \item{chipType}{@character, chip type.}
#   \item{...}{arguments passed to 
#      \code{AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(...)}.}
#   \item{rootPath}{@character, root path for output files.}
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
# \note{
#   The implementation requires \pkg{aroma.affymetrix} only for the
#   inference of a @see "aroma.core::UnitNamesFile".  However it is tied
#   to Affymetrix data anyway because Birdseed only supports Affymetrix data.
# }
#
#*/###########################################################################
setMethodS3("exportGenotypeCallsAndConfidenceScores", "BroadBirdseedGenotypeTcgaDataFile", function(this, dataset, chipType, ..., rootPath="callData", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet, length=c(1,1));

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  path <- file.path(rootPath, dataSet, chipType);

  fullname <- getFullName(this);
  verbose && enter(verbose, "Exporting genotype calls and confidence scores for sample ", fullname);
  # Nothing to do ?
  genoFilename <- sprintf("%s,genotypes.acf", fullname);
  genoPathname <- Arguments$getWritablePathname(genoFilename, path=path, mustNotExist=FALSE);
  confFilename <- sprintf("%s,confidenceScores.asb", fullname);
  confPathname <- Arguments$getWritablePathname(confFilename, path=path, mustNotExist=FALSE);
  
  if (isFile(genoPathname) && isFile(confPathname)) {
    acf <- AromaUnitSignalBinaryFile(genoPathname);
    verbose && cat(verbose, "Already exported: ", getPathname(acf));
    asf <- AromaUnitSignalBinaryFile(confPathname);
    verbose && cat(verbose, "Already exported: ", getPathname(asf));
  } else {
    require("aroma.affymetrix") || 
      throw("Package not loaded: aroma.affymetrix");
    # Affy-specific
    verbose && enter(verbose, "Retrieving UnitNamesFile from chipType");
    unf <- AffymetrixCdfFile$byChipType(chipType, verbose=verbose);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Defining file header");
    srcFile <- list(filename = fullname,
                    fileSize = getFileSize(this),
                    checkSum = getChecksum(this));
    verbose && exit(verbose);
      
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Transforming data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Extracting calls and confidence scores");
    calls <- extractCalls(this);
    verbose && str(verbose, calls);
    verbose && print(verbose, summary.factor(calls));
    
    conf <- 1-extractConfidenceScores(this); # input data are p-values...
    verbose && str(verbose, conf);
    
    verbose && enter(verbose, "Ordering units according to UnitNamesFile");
    data <- readDataFrame(df);
    units <- indexOf(unf, names=data[, "CompositeElement REF"]);
    stopifnot(!anyMissing(units));
    verbose && exit(verbose);
    
    rm(data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Extracting genotype calls");
    filename <- sprintf("%s,genotypes.acf", fullname);
    
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE);
    
    if (isFile(pathname)) {
      acf <- AromaUnitSignalBinaryFile(pathname);
      verbose && cat(verbose, "Already exported: ", getPathname(acf));
    } else {
      
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Writing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pathnameT <- sprintf("%s.tmp", pathname);
      pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);
      on.exit({
        file.remove(pathnameT);
      }, add=TRUE);
      
      acf <- AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(unf, filename=pathnameT, overwrite=FALSE, ...);
      
      updateGenotypes(acf, units=units, calls=calls, encoding="birdseed", verbose=verbose);
      verbose && exit(verbose);
      
      footer <- readFooter(acf);
      footer$srcFile <- srcFile;
      writeFooter(acf, footer);
      
      # Rename temporary file
      res <- file.rename(pathnameT, pathname);
      if (!isFile(pathname)) {
        throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
      }
      if (isFile(pathnameT)) {
        throw("Failed to remove temporary file: ", pathnameT);
      }
    }
    
    verbose && enter(verbose, "Extracting confidence scores");
    filename <- sprintf("%s,confidenceScores.asb", fullname);
    
    srcFile$filename <- filename;
    
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE);
    
    if (isFile(pathname)) {
      asf <- AromaUnitSignalBinaryFile(pathname);
      verbose && cat(verbose, "Already exported: ", getPathname(asf));
    } else {
      pathnameT <- sprintf("%s.tmp", pathname);
      pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);
      on.exit({
        file.remove(pathnameT);
      }, add=TRUE);
      
      asf <- AromaUnitSignalBinaryFile$allocateFromUnitNamesFile(unf, filename=pathnameT, types="double", size=4, signed=TRUE, overwrite=FALSE, ...);
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
      res <- file.rename(pathnameT, pathname);
      if (!isFile(pathname)) {
        throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
      }
      if (isFile(pathnameT)) {
        throw("Failed to remove temporary file: ", pathnameT);
      }
    }
    verbose && exit(verbose);
  }
  
  verbose && exit(verbose);

  invisible(list(acf, asf));
})


############################################################################
# HISTORY:
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
