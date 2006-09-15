setMethodS3("print", "GenericSummary", function(x, ..., collapse="\n") {
  # To please R CMD check
  this <- x;

  s <- paste(this, collapse=collapse);
  cat(s, collapse, sep="");
})


############################################################################
# HISTORY:
# 2006-08-11
# o Created.
############################################################################
