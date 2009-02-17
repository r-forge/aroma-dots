###########################################################################/**
# @RdocClass AromaUnitTotalCnBinaryFile
#
# @title "The AromaUnitTotalCnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitTotalCnBinaryFile is a @see "AromaUnitSignalBinaryFile".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitSignalBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/########################################################################### 
setConstructorS3("AromaUnitTotalCnBinaryFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), "AromaUnitTotalCnBinaryFile"
  );
})


setMethodS3("extractRawCopyNumbers", "AromaUnitTotalCnBinaryFile", function(this, chromosome, range=NULL, units=NULL, ..., verbose=FALSE) {
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(this)));
    units <- sort(unique(units));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Extracting RawCopyNumbers");

  name <- getFullName(this);
  verbose && cat(verbose, "Name: ", name);

  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && enter(verbose, "Identifying units on chromosome");
  ugp <- getAromaUgpFile(this, ..., verbose=less(verbose,50));
  verbose && print(verbose, ugp);
  units2 <- getUnitsAt(ugp, chromosome=chromosome, range=range, ..., 
                                            verbose=less(verbose,5));
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units2);

  # Keeping only a subset of units?
  if (!is.null(units)) {
    verbose && enter(verbose, "Keeping only units of interest");
    keep <- is.element(units2, units);
    verbose && cat(verbose, "Keeping:");
    verbose && summary(verbose, keep);
    units2 <- units2[keep];
    rm(keep);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    verbose && exit(verbose);
  }
  units <- units2;
  rm(units2);


  verbose && cat(verbose, "Genomic positions:");
  pos <- getPositions(ugp, units=units);
  verbose && str(verbose, pos);  
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting data");
  M <- extractMatrix(this, units=units, drop=TRUE, verbose=less(verbose,5));
  verbose && str(verbose, M);
  rawCNs <- RawCopyNumbers(x=pos, cn=M, chromosome=chromosome, name=name);

  # Add annotation data
  rawCNs$platform <- getPlatform(this);
  rawCNs$chipType <- getChipType(this);
  rawCNs$fullname <- getFullName(this);

  verbose && exit(verbose);

  verbose && exit(verbose);

  rawCNs;
})



############################################################################
# HISTORY:
# 2009-02-17
# o Added argument 'units' to extractRawCopyNumbers().
# 2009-02-16
# o Now extractRawCopyNumbers() also includes the full (sample) name.
# 2008-06-12
# o Now extractRawCopyNumbers() adds annotation data to the returned object,
#   i.e. platform, chipType, and fullname.
# 2008-05-21
# o Added extractRawCopyNumbers().
# 2008-05-11
# o Created.
############################################################################
