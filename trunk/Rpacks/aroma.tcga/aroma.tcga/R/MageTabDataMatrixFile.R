setConstructorS3("MageTabDataMatrixFile", function(..., .verify=TRUE) {
  this <- extend(MageTabFile(..., .verify=FALSE), "MageTabDataMatrixFile");

#  if (.verify)
#    verify(this, ..., verbose=verbose);
  this;
})



setMethodS3("getHeader", "MageTabDataMatrixFile", function(this, ..., force=FALSE) {
  hdr <- this$.fileHeader;
  if (force || is.null(hdr)) {
    hdr <- readRawHeader(this, ...);
    if (hasColumnHeader(this)) {
      header <- list();
      for (kk in 1:2) {
        tmp <- hdr$topRows[[kk]];
        header[[tmp[1]]] <- tmp[-1];
      }
      hdr$dataMatrixHeader <- header;
      colnames <- sprintf("%s,%s", header[[1]], header[[2]]);
      colnames <- gsub("(^,|,$)", "", colnames);
      colnames <- c(names(header)[2], colnames);
      hdr$columns <- colnames;
      hdr$skip <- hdr$skip + 1L;
    }
    this$.fileHeader <- hdr;
  }
  hdr;
})


setMethodS3("getReadArguments", "MageTabDataMatrixFile", function(this, ..., colClassPatterns=c("*"="character", "(Signal)$"="double")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE);



############################################################################
# HISTORY:
# 2009-05-18
# o Added to aroma.tcga.
# 2009-04-19
# o Fixed the column names.
# 2009-04-18
# o Created.
############################################################################
