###########################################################################/**
# @RdocClass AromaTotalCnBinaryFile
#
# @title "The AromaTotalCnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaTotalCnBinaryFile is a @see "AromaTabularBinaryFile".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/########################################################################### 
setConstructorS3("AromaTotalCnBinaryFile", function(...) {
  extend(AromaSignalBinaryFile(...), "AromaTotalCnBinaryFile"
  );
})


setMethodS3("extractRawCopyNumbers", "AromaTotalCnBinaryFile", function(this, chromosome, range=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Extracting RawCopyNumbers");

  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && enter(verbose, "Identifying units on chromosome");
  ugp <- getAromaUgpFile(this, ..., verbose=less(verbose,50));
  verbose && print(verbose, ugp);
  units <- getUnitsAt(ugp, chromosome=chromosome, range=range, ..., verbose=less(verbose,5));
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);  

  verbose && cat(verbose, "Genomic positions:");
  pos <- getPositions(ugp, units=units);
  verbose && str(verbose, pos);  
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting data");
  M <- extractMatrix(this, rows=units, drop=TRUE, verbose=less(verbose,5));
  verbose && str(verbose, M);  
  rawCNs <- RawCopyNumbers(x=pos, cn=M, chromosome=chromosome);
  verbose && exit(verbose);

  verbose && exit(verbose);

  rawCNs;
})



############################################################################
# HISTORY:
# 2008-05-21
# o Added extractRawCopyNumbers().
# 2008-05-11
# o Created.
############################################################################
