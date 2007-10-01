setMethodS3("drawExtraAnnotations", "default", function(fit, ...) {
  rawCNs <- extractRawCopyNumbers(fit);
  rawCNs <- getCNs(rawCNs);

  sd <- sd(rawCNs, na.rm=TRUE);
  mad <- mad(rawCNs, na.rm=TRUE);

  text <- substitute(paste(hat(sigma)==sd, ", ", hat(sigma)["MAD"]==mad, sep=""), 
                   list(sd=sprintf("%.3g", sd), mad=sprintf("%.3g", mad)));
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
# 2007-09-29
# o Added default annotations with sd and MAD estimates.
# 2007-09-04
# o Created.
############################################################################
