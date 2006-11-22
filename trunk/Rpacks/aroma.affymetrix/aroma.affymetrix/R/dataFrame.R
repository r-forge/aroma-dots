setMethodS3("dataFrame", "default", function(colClasses, nrow=1, ...) {
  df <- vector("list", length(colClasses));
  names(df) <- names(colClasses);
  for (kk in seq(along=df)) {
    df[[kk]] <- vector(colClasses[kk], nrow);
  }

  attr(df, "row.names") <- seq(length=nrow);
  class(df) <- "data.frame";
  df;
})

##############################################################################
# HISTORY:
# 2006-11-22
# o Created.
##############################################################################
