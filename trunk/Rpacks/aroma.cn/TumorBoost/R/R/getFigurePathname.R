getFigurePathname <- function(name, tags=NULL, ext="png", path="figures", ...) {
  fullname <- paste(c(name, tags), collapse=",");
  filename <- sprintf("%s.%s", fullname, ext);
  pathname <- Arguments$getWritablePathname(filename, path=path);
  pathname;
} # getFigurePathname()


############################################################################
# HISTORY:
# 2009-03-02 [HB]
# o Created.
############################################################################
