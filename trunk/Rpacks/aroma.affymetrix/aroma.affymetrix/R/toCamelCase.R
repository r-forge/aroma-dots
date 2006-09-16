setMethodS3("toCamelCase", "default", function(s, capitalize=FALSE, split="[ \t]+", ...) {
  s <- strsplit(s, split=split);
  s <- lapply(s, FUN=function(s) {
    paste(capitalize(tolower(s)), collapse="")
  });
  s <- unlist(s);
  if (!capitalize)
    s <- decapitalize(s);
  s;
})

############################################################################
# HISTORY:
# 2006-09-15
# o Created.  Will probably end up in R.utils some day.
############################################################################
