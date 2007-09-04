setMethodS3("drawExtraAnnotations", "default", function(fit, ...) {
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
# 2007-09-04
# o Created.
############################################################################
