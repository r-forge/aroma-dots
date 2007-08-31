.patchRowMedians <- function(...) {
  # 1) rowMedians() of Biobase v1.15.27 and above got it all.
  if (require("Biobase")) {
    ver <- packageDescription("Biobase")$Version;
    if (compareVersion(ver, "1.15.27") >= 0) {
      # Nothing to do.
      return();
    }
  }

  # 2) rowMedians() of R.native v0.1.2 and above got it all.
  if (require("R.native")) {
    ver <- packageDescription("R.native")$Version;
    if (compareVersion(ver, "0.1.2") >= 0) {
      # Nothing to do.
      return();
    }
  }

  # As a final resort, we use apply():
  setMethodS3("rowMedians", "matrix", function(X, na.rm=FALSE, ...) {
    # About 3-10 times slower than Biobase::rowMedians()
    base::apply(X, MARGIN=1, FUN=median, na.rm=na.rm, ...);
  })
} # .patchRowMedians()


############################################################################
# HISTORY:
# 2007-08-30
# o Added to get a consistent rowMedians() function across different 
#   versions of R, Biobase, and with/or without R.native.  Eventually,
#   all that will be need is Biobase::rowMedians(), but until then we
#   have to wait.
############################################################################
