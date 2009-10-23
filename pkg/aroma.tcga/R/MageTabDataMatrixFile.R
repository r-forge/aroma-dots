setConstructorS3("MageTabDataMatrixFile", function(..., .verify=TRUE) {
  this <- extend(MageTabFile(..., .verify=FALSE), "MageTabDataMatrixFile");

#  if (.verify)
#    verify(this, ..., verbose=verbose);
  this;
})



setMethodS3("getHeader", "MageTabDataMatrixFile", function(this, ..., strict=FALSE, force=FALSE) {
  hdr <- this$.fileHeader;
  if (force || is.null(hdr)) {
    hdr <- readRawHeader(this, ...);
    if (hasColumnHeader(this)) {
      header <- list();
      for (kk in 1:2) {
        tmp <- hdr$topRows[[kk]];
        header[[tmp[1]]] <- tmp[-1];
      }

      ns <- sapply(header, FUN=length);
      if (ns[1] < ns[2]) {
        if (strict) {
          throw("Broken header: ns[1] < ns[2]");
        }
        # Very ad hoc; solves a specific case with the Broad. /HB 2009-09-25
        dn <- ns[2] - ns[1];
        header[[1]] <- c(rep("", length.out=dn), header[[1]]);
      } else if (ns[1] > ns[2]) {
        throw("Broken header: ns[1] < ns[2]");
      }

      hdr$dataMatrixHeader <- header;
      colnames <- sprintf("%s,%s", header[[1]], header[[2]]);
      colnames <- gsub("(^,|,$)", "", colnames);
      colnames <- gsub("(Rtum/Rnorm)", "Ratio", colnames, fixed=TRUE);
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
# 2009-09-24
# o Added argument 'strict=FALSE' to getHeader() of MageTabDataMatrixFile.
#   If the number of 'Hybridization REF' and 'CompositeElement REF' differ,
#   which is not a valid MAGE_TAB format, and strict=FALSE, then the
#   'Hybridization REF' is padded with "" from the beginning.
# 2009-08-20
# o Now getHeader() of MageTabDataMatrixFile also recognizes columns named
#   '(Rtum/Rnorm)'.
# 2009-05-18
# o Added to aroma.tcga.
# 2009-04-19
# o Fixed the column names.
# 2009-04-18
# o Created.
############################################################################
