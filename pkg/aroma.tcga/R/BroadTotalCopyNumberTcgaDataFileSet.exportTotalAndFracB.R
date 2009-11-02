setMethodS3("exportTotalAndFracB", "BroadTotalCopyNumberTcgaDataFileSet", function(this, unf, ..., rootPath="totalAndFracBData", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  dataSet <- getFullName(this);
  verbose && cat(verbose, "Full data set name: ", dataSet);

  chipType <- getChipType(unf, fullname=FALSE);
  
  path <- file.path(rootPath, dataSet, chipType);
  path <- Arguments$getWritablePath(path);


  for (kk in seq(this)) {
    df <- getFile(this, kk);
    
    verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                           kk, getName(df), length(this)));

    dsKK <- exportTotalAndFracB(df, dataSet=dataSet, unf=unf, 
                           rootPath=rootPath, ..., verbose=less(verbose, 5));
    verbose && cat(verbose, "Exported data set:");
    verbose && print(verbose, dsKK);

    verbose && exit(verbose);
  } # for (cc ...)

  # Validating everything
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
