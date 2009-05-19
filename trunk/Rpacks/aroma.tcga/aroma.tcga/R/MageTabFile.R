setConstructorS3("MageTabFile", function(..., .verify=TRUE) {
  this <- extend(TabularTextFile(..., .verify=FALSE), "MageTabFile");

  if (.verify)
    verify(this, ..., verbose=verbose);
  this;
})



############################################################################
# HISTORY:
# 2009-05-18
# o Added to aroma.tcga.
# 2009-04-19
# o Fixed the column names.
# 2009-04-18
# o Created.
############################################################################
