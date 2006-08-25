setMethodS3("print", "GenericSummary", function(object, ..., collapse="\n") {
  s <- paste(object, collapse=collapse);
  cat(s, collapse, sep="");
})


############################################################################
# HISTORY:
# 2006-08-11
# o Created.
############################################################################
