###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the RMA 
# chip-effect estimates as estimated by affyPLM.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2007-06-20
# Last modified: 2008-02-19
#
# Data set:
#  rawData/
#   Affymetrix-HeartBrain/
#    HG-U133_Plus_2/
#     u1332plus_ivt_cerebellum_A.CEL [13555904 bytes]
#     u1332plus_ivt_cerebellum_B.CEL [13550687 bytes]
#     u1332plus_ivt_cerebellum_C.CEL [13551860 bytes]
#     u1332plus_ivt_heart_A.CEL      [13554731 bytes]
#     u1332plus_ivt_heart_B.CEL      [13553255 bytes]
#     u1332plus_ivt_heart_C.CEL      [13551203 bytes]
#  Source: Affymetrix Tissue samples, 2007.  http://www.affymetrix.com/
#  support/technical/sample_data/hugene_1_0_array_data.affx
###########################################################################

library(aroma.affymetrix);
library(affyPLM);          # fitPLM()

verbose <- Arguments$getVerbose(-8);
timestampOn(verbose);

doPlot <- TRUE;
saveImg <- TRUE;

# ----------------------------------
# RMA estimates by aroma.affymetrix 
# ----------------------------------
verbose && enter(verbose, "RMA by aroma.affymetrix");

cs <- AffymetrixCelSet$byName("Affymetrix-HeartBrain", 
                              chipType="HG-U133_Plus_2");

# RMA background correction
bc <- RmaBackgroundCorrection(cs);
csBC <- process(bc, verbose=verbose);

# RMA quantile normalization
qn <- QuantileNormalization(csBC, typesToUpdate="pm");
csN <- process(qn, verbose=verbose);

# RMA probe summarization
plm <- RmaPlm(csN);
fit(plm, verbose=verbose);

# Extract chip effects on the log2 scale
ces <- getChipEffectSet(plm);
theta <- extractMatrix(ces, returnUgcMap=TRUE);
theta <- log2(theta);
ugcMap <- attr(theta, "unitGroupCellMap");
rownames(theta) <- getUnitNames(getCdf(ces), ugcMap[,"unit"]);

verbose && exit(verbose);


# ------------------------
# RMA estimate by affyPLM
# ------------------------
verbose && enter(verbose, "RMA by aroma.affymetrix");

raw <- ReadAffy(filenames=getPathnames(cs));
fit <- fitPLM(raw, verbos=9);
theta0 <- coefs(fit);

verbose && exit(verbose);


# --------------------------------
# Compare the two implementations
# --------------------------------
# Reorder the aroma.affymetrix estimates
o <- match(rownames(theta0), rownames(theta));
theta <- theta[o,];

e <- (theta - theta0);
stopifnot(mean(as.vector(e^2)) < 1e-4);
stopifnot(sd(as.vector(e^2)) < 1e-4);

if (doPlot) {
  if (saveImg) {
    pngDev <- findPngDevice();
    pngDev("replication-affyPLM,fitPLM.png", width=640, height=640);
  }

  layout(matrix(1:9, ncol=3, byrow=TRUE));

  xlab <- expression(log[2](theta[affyPLM]));
  ylab <- expression(log[2](theta[aroma.affymetrix]));
  for (kk in seq(length=ncol(theta))) {
    main <- colnames(theta)[kk];
    plot(theta0[,kk], theta[,kk], pch=".", xlab=xlab, ylab=ylab, main=main);
    abline(0,1, col="blue");
  }

  xlab <- expression(log[2](theta[aroma.affymetrix]/theta[affyPLM]));
  plotDensity(e, xlab=xlab);

  if (saveImg)
    dev.off();
}
