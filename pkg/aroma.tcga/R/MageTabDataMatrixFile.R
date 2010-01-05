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

      # Assume two header rows
      # Example:
      # Hybridization REF<tab>TCGA-12-0620-01A-01D-0310-01
      # CompositeElement REF<tab>Call<tab>Confidence
      header <- hdr$topRows[1:2];

      # However, some old TCGA data files only having one header row, i.e.
      # they are plain old tab-delimited text files.
      # Try to detect such old ones by looking for columns that are numbers.
      hasNums <- sapply(header, FUN=function(x) {
        suppressWarnings({
          y <- as.double(x);
        });
        isNum <- (is.finite(y) & (x == y));
        any(isNum);
      });
      # If 2nd header row has pure numbers, assume that is not a header row.
      if (hasNums[2]) {
        if (strict) {
          throw("Not a valid MAGE-TAB header (2nd row contains numbers): ", getPathname(this));
        }
        header[[2]] <- NULL;
        isPlainCSV <- TRUE;
      } else {
        isPlainCSV <- FALSE;
      }

      # Sanity check; don't allow numbers in the first header row
      if (hasNums[1]) {
        throw(sprintf("File format error of %s ('%s'). First header row contains columns that are pure numbers: %s", 
                            class(this), getPathname(this), paste(header[[1]], collapse=", ")));
      }

      # Is it a standard tab-delimited text file?
      if (isPlainCSV) {
        # (a) A standard tab-delimited text file
        colnames <- header[[1]];
      } else {
        # (b) A MAGE TAB file

        # Sanity checks
        ns <- sapply(header, FUN=length);
        if (ns[1] < ns[2]) {
          if (strict) {
            throw("Not a valid MAGE-TAB header (ns[1] < ns[2]): ", getPathname(this));
          }
          # Very ad hoc; solves a specific case with the Broad. /HB 2009-09-25
          dn <- ns[2] - ns[1];
          header[[1]] <- c(rep("", length.out=dn), header[[1]]);
        } else if (ns[1] > ns[2]) {
          throw("Not a valid MAGE-TAB header (ns[1] > ns[2]): ", getPathname(this));
        }

        # Drop 'Hybridization REF'/'CompositeElement REF' entries
        names <- sapply(header, FUN=function(s) s[1]);
        header <- lapply(header, FUN=function(s) s[-1]);
        names(header) <- names;

        colnames <- sprintf("%s,%s", header[[1]], header[[2]]);
        colnames <- gsub("(^,|,$)", "", colnames);
        colnames <- c(names(header)[2], colnames);

        hdr$skip <- hdr$skip + 1L;

        # AD HOC/FIX.  /HB 2009-08-20
        colnames <- gsub("(Rtum/Rnorm)", "Ratio", colnames, fixed=TRUE);
      }

      # Paste headers      
      hdr$dataMatrixHeader <- header;
      hdr$columns <- colnames;
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
# 2010-01-05
# o ROBUSTNESS: Some of the TCGA data files are not MAGE-TAB files, but
#   rather regular tab-delimited files.  Updated getHeader() of 
#   MageTabDataMatrixFile to detect this.
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
