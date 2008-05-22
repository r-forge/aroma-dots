setMethodS3("extractRawCopyNumbers", "profileCGH", function(object, ...) {
  pv <- object$profileValues;
  chromosome <- unique(pv$Chromosome);
  RawCopyNumbers(cn=pv$LogRatio, x=pv$PosBase, chromosome=chromosome);
})


setMethodS3("extractRawCopyNumbers", "DNAcopy", function(object, ...) {
  data <- object$data;
  chromosome <- unique(data$chrom);
  RawCopyNumbers(cn=data[[3]], x=data$maploc, chromosome=chromosome);
})




############################################################################
# HISTORY:
# 2008-05-21
# o Now extractRawCopyNumbers() adds 'chromosome' to the returned object.
# 2008-05-17
# o Extracted RawCopyNumbers.FOREIGN.R after moving RawCopyNumbers.R
#   to aroma.core.  
# 2008-03-31
# o Put recently added sd() and mad() into estimateStandardDeviation().
# 2008-03-10
# o Added standard deviation estimator sd() and mad() which my default
#   uses a first-order difference variance estimator.
# 2007-08-22
# o Created.  Need a generic container for holding copy number data and
#   to plot them nicely.
############################################################################
