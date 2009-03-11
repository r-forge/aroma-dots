rspCapture <- function(..., wrapAt=80, collapse="\n", FUN=print) {
  x <- capture.output(FUN(...));

  # Wrap long lines?
  if (!is.null(wrapAt)) {
    nok <- (nchar(x) > wrapAt);
    if (any(nok)) {
      x <- as.list(x);
      x[nok] <- lapply(x[nok], FUN=function(s) {
        res <- NULL;
        while(nchar(s) > 0) {
          res <- c(res, substr(s, 1, wrapAt));
          s <- substr(s, wrapAt+1, nchar(s));
        }
        res;
      });
      x <- unlist(x);
    }
  }
  
  # Concatenate rows?
  if (!is.null(collapse)) {
    x <- paste(x, collapse=collapse);
  }

  x;
} # rspCapture()


###############################################################################
# HISTORY:
# 2009-03-11
# o Added argument FUN=print.
# 2009-02-25
# o Created.
###############################################################################
