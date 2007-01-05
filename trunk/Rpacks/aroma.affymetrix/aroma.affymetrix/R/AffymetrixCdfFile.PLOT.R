setMethodS3("stextChipType", "AffymetrixCdfFile", function(this, side=4, fmtstr="%s", pos=1, cex=0.7, col="darkgray", ...) {
  stext(side=side, text=sprintf(fmtstr, getChipType(this)), pos=pos, cex=cex, col=col, ...);
}, private=TRUE)



############################################################################
# HISTORY:
# 2006-09-16
# o Added getGenomeInformation() and stextChipType().
############################################################################
