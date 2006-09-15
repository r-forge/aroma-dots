setMethodS3("normalizeGenotypeCone", "matrix", function(y, avg=median, targetAvg=2200, ...) {
  fit <- fitGenotypeCone(y, ...);
  x <- backtransformGenotypeCone(y, fit=fit);
#  x <- normalizeAverage(x, avg=avg, targetAvg=targetAvg);
  x;
}) # normalizeGenotypeCone()



############################################################################
# HISTORY:
# 2006-05-08
# o Created.
############################################################################
