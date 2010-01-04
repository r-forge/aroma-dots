###########################################################################/**
# @set "class=HudsonAlphaGenotypeTcgaDataFile"
# @RdocMethod exportGenotypeCalls
#
# @title "Exports genotype calls and confidence scores"
#
#
# \description{
#    Export genotype calls and corresponding confidence scores for one sample
#    from a BeadStudio data file to an 
#     @see "aroma.core::AromaUnitGenotypeCallFile" (genotype calls).
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
#   A @see "aroma.core::UnitNamesFile" is inferred from "chipType".
#   Genotyping units (SNPs) are ordered according to this 
#   @see "aroma.core::UnitNamesFile".
# }
#
#*/###########################################################################
setMethodS3("exportGenotypeCalls", "HudsonAlphaGenotypeTcgaDataFile", function(this, dataSet, unf, ..., rootPath="callData", force=FALSE, verbose=FALSE) {
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
  tags <- c("genotypes");

  # Unit indices (to be inferred)
  units <- NULL;

  # Export total signals
  allColumnNames <- getColumnNames(this);
  pattern <- ",genotype$";
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

    filename <- sprintf("%s.acf", fullname);
    pathname <- file.path(path, filename);
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
    data <- extractCalls(this, sampleNames=sampleName, drop=FALSE, verbose=less(verbose,10));
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
    df <- AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(unf, 
             filename=pathnameT, path=NULL, overwrite=force, ...);
      
    footer <- readFooter(df);
    footer$srcFile <- srcFile;
    writeFooter(df, footer);
    verbose && exit(verbose);

    verbose && enter(verbose, "Write signals");
    updateGenotypes(df, units=units, calls=data, encoding="generic", 
                    verbose=verbose);
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
    df <- AromaUnitGenotypeCallFile(pathname);
    verbose && cat(verbose, "Exported data file:");
    verbose && print(verbose, df);

    verbose && exit(verbose);
  } # for (cc ...)

  ds <- AromaUnitGenotypeCallSet$byPath(path);
  verbose && print(verbose, ds);

  verbose && exit(verbose);

  invisible(ds);
})


############################################################################
# HISTORY:
# 2009-12-05
# o Created.
############################################################################
