setConstructorS3("AffymetrixNetAffxCsvFile", function(..., .verify=TRUE) {
  this <- extend(AffymetrixCsvFile(..., .verify=FALSE), "AffymetrixNetAffxCsvFile");

  if (.verify)
    verify(this, ...);
  this;
})

setMethodS3("findByChipType", "AffymetrixNetAffxCsvFile", function(static, chipType, tags=".*", pattern=sprintf("^%s%s([.]|_)annot[.](csv|CSV)$", chipType, tags), ...) {
  findByChipType.AffymetrixCsvFile(static, chipType=chipType, pattern=pattern, ...);
}, static=TRUE, protected=TRUE)


setMethodS3("readDataUnitChromosomePosition", "AffymetrixNetAffxCsvFile", function(this, colClassPatterns=c("*"="NULL", "^probeSetID$"="character", "^chromosome$"="character", "^(physicalPosition|chromosomeStart)$"="character"), con=NULL, ...) {
  data <- readData(this, colClassPatterns=colClassPatterns, camelCaseNames=TRUE, ...);

  # Convert chromosome strings to integers
  cc <- grep("^chr", colnames(data))[1];
  map <- c(X=23, Y=24, Z=25);
  for (kk in seq(along=map)) {
    data[[cc]] <- gsub(names(map)[kk], map[kk], data[[cc]]);
  }
  data[[cc]] <- as.integer(data[[cc]]);
  gc <- gc();

  # Convert positions to integers
  cc <- grep("(p|P)os", colnames(data));
  if (length(cc) == 0)
    cc <- grep("(s|S)tart", colnames(data));
  data[[cc]] <- as.integer(data[[cc]]);
  gc <- gc();

  colnames(data) <- c("unitName", "chromosome", "position");
  attr(data, "header") <- NULL;

  data;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
