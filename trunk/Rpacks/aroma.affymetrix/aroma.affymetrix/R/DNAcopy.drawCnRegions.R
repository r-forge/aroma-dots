setMethodS3("drawCnRegions", "DNAcopy", function(this, ...) {
  output <- this$output;
  starts <- output[["loc.start"]];
  ends <- output[["loc.end"]];
  means <- output[["seg.mean"]];
  nbrOfSegments <- length(means);
  xx <- cbind(starts, ends);
  yy <- cbind(means, means);
  for (rr in 1:nrow(xx)) {
    lines(x=xx[rr,], y=yy[rr,], ...);
  };
})

##############################################################################
# HISTORY:
# 2007-08-20
# o Created.
##############################################################################
