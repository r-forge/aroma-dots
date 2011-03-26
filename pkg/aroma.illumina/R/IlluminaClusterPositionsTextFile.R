setConstructorS3("IlluminaClusterPositionsTextFile", function(...) {
  extend(TabularTextFile(...), "IlluminaClusterPositionsTextFile");
})


setMethodS3("readCentroidMatrix", "IlluminaClusterPositionsTextFile", function(this, ..., what=c("T", "R")) {
  # Argument 'what':
  what <- match.arg(what);

  fieldPattern <- sprintf("(AA|AB|BB)_%s_Mean", what);
  patterns <- c(reporterId="character", "(AA|AB|BB)_._Mean"="numeric");
  names(patterns)[2] <- fieldPattern;

  data <- readDataFrame(this, colClassPatterns=patterns, ...);
  mu <- as.matrix(data[,-1,drop=FALSE]);
  rownames(mu) <- data[,1,drop=TRUE];
  colnames(mu) <- gsub(fieldPattern, "\\1", colnames(mu));

  # Sanity check
  stopifnot(all(colnames(mu) == c("AA", "AB", "BB")));

  mu;
}) # readCentroidMatrix()


############################################################################
# HISTORY:
# 2011-03-25
# o Added readCentroidMatrix().
# o Created.  Just a stub.
############################################################################
