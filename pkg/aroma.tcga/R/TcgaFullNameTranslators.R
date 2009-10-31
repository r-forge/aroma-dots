setConstructorS3("TcgaFullNameTranslators", function(...) {
  extend(Object(...), "TcgaFullNameTranslators");
})

setMethodS3("getBarcodeIdTranslators", "TcgaFullNameTranslators", function(static, ...) {
  require("gsubfn") || throw("Package not loaded: gsubfn");

  # Example: 
  # TCGA-16-0849-01A-01D-0384-01 => TCGA-16-0849,01A,01D-0384-01
  pattern <- BiospecimenCoreResource$getBarcodePattern("aliqoutBarcode");
  fnt <- function(names, ...) {
    # Workaround: strapply() drops elements not matching the pattern
    idxs <- grep(pattern, names);
    names[idxs] <- strapply(names[idxs], pattern=pattern, FUN=function(...) {
      x <- list(...);
      sprintf("%s,%s,%s-%s", x[2], x[6], x[9], x[12]);
    });
    unlist(names, use.names=FALSE);
  } # fnt()

  res <- list(
    "sampleName,tumorType,rest" = fnt    
  );

  res;
}, static=TRUE)


setMethodS3("getBarcodeIdTranslator", "TcgaFullNameTranslators", function(static, name, ...) {
  fnList <- getBarcodeIdTranslators(static, ...);
  fnList[[name]];
}, static=TRUE)



############################################################################
# HISTORY:
# 2009-10-30
# o BUG FIX: getBarcodeIdTranslators() of TcgaFullNameTranslators would
#   drop fullnames that did not match the pattern.
# 2009-10-22
# o Created.
############################################################################
