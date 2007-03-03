setConstructorS3("AromaCgpFile", function(...) {
  extend(AromaGenomePositionFile(...), "AromaCgpFile")
})

setMethodS3("getFilenameExtension", "AromaGenomePositionFile", function(static, ...) {
  "cgp";
}, static=TRUE)

setMethodS3("nbrOfCells", "AromaCgpFile", function(this, ...) {
  nbrOfElements(this, ...);
})


setMethodS3("createFromCdf", "AromaCgpFile", function(static, cdf, ...) {
  chipType <- getChipType(cdf);
  create(static, chipType=chipType, nbrOfElements=nbrOfCells(cdf), ...);
}, static=TRUE)



setMethodS3("getCellsAt", "AromaCgpFile", function(this, ...) {
  getElementsAt(this, ...);
})

setMethodS3("importFromAffymetrixCsvFile", "AromaCgpFile", function(this, csv, ...) {
  # To do
})


############################################################################
# HISTORY:
# 2007-03-03
# o Created from AromaUgpFile.R.
############################################################################
