setMethodS3("drawCnRegions", "profileCGH", function(this, xscale=1, ...) {
  # Get data
  pv <- this$profileValues;

  # Order data along chromosome(s)
  o <- order(pv$Chromosome, pv$PosBase);
  pv <- pv[o, ];

  # Number of data points
  n <- length(pv[, 1]);

  posMax <- max(pv$PosBase) + 1;
  pos <- pv$PosBase[1:(n-1)];
  posNext <- pv$PosBase[2:n];

  interPos <- pos + (posNext - pos)/2;
  interPos <- c(0, interPos, posMax);
  smtStart <- pv[, "Smoothing"][1];
  smtEnd <- pv[, "Smoothing"][n];
  smt1 <- pv[, "Smoothing"][1:(n-1)];
  smt1 <- c(smtStart, smt1, smtEnd);
  smt2 <- pv[, "Smoothing"][2:n];
  smt2 <- c(smtStart, smt2, smtEnd);
  datasmt <- data.frame(posBase=c(interPos, interPos), smoothing=c(smt1, smt2));
  datasmt <- unique(datasmt);
  datasmt <- datasmt[order(datasmt$posBase), ];
#  posBase <- xscale * datasmt$posBase;

  lines(datasmt$smoothing ~ datasmt$posBase, ...);
})


##############################################################################
# HISTORY:
# 2007-08-20
# o Created.
##############################################################################
