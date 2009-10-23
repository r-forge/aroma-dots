setConstructorS3("MscnExport", function(...) {
  extend(Object(), "MscnExport");
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fullnames translator for data files
# Pattern: <samplePattern>,ref=<samplePattern>,log2ratio,total
#
# Example:
# Turn this fullname:
# TCGA-25-1323,01A,01D-0452-01,ref=TCGA-25-1323,10A,01D-0452-01,log2ratio,total
# into this one:
# TCGA-25-1323,01Avs10A,01D-0452-01vs01D-0452-01,log2ratio,total
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getFullNameTranslator", "MscnExport", function(static, ...) {
  require("gsubfn") || throw("Package not loaded: gsubfn");

  patterns <- BiospecimenCoreResource$getBarcodePatterns();
  samplePattern <- with(patterns, sprintf("(%s),(%s),(%s-%s)", 
                        patient, sampleId, portionId, plateBarcode));
  pattern <- sprintf("^%s,ref[=]%s,(.*)", samplePattern, samplePattern);

  translator <- function(names, ...) {
    names <- strapply(names, pattern=pattern, FUN=function(...) {
      x <- unlist(list(...));
      n <- length(x);
      sprintf("%s,%svs%s,%svs%s,%s", x[1], x[5], x[17], x[8], x[20], x[n]);
    });
    names;
  } # translator()

  translator;
}, static=TRUE)


############################################################################
# HISTORY:
# 2009-10-03
# o Created.
############################################################################
