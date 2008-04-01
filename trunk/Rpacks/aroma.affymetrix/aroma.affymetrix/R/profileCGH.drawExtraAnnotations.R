setMethodS3("drawExtraAnnotations", "default", function(fit, ...) {
  rawCNs <- extractRawCopyNumbers(fit);
  rawCNs <- getCNs(rawCNs);
  sd <- estimateStandardDeviation(rawCNs);
  text <- substitute(hat(sigma)==sd, list(sd=sprintf("%.3g", sd)));
  stext(text=text, side=3, pos=0.5, line=-2); 
}, protected=TRUE);


setMethodS3("drawExtraAnnotations", "profileCGH", function(fit, ...) {
  sdEst <- fit$SigmaC$Value;
  if (!is.null(sdEst)) {
    text <- substitute(hat(sigma)==x, list(x=sprintf("%.3g", sdEst)));
    stext(text=text, side=3, pos=0.5, line=-2);
  }
}, protected=TRUE);


############################################################################
# HISTORY:
# 2008-03-31
# o Now the standard deviation across all CNs in a chromosome is calculated
#   using a robust first-order difference estimator, which will make the 
#   estimate much less affected by copy-number changes.
# 2007-09-29
# o Added default annotations with sd and MAD estimates.
# 2007-09-04
# o Created.
############################################################################
