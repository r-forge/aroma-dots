setConstructorS3("BroadBirdseedGenotypeTcgaDataFile", function(...) {
  this <- extend(TcgaDataFile(...), "BroadBirdseedGenotypeTcgaDataFile");
  this;
})


setMethodS3("getExtensionPattern", "BroadBirdseedGenotypeTcgaDataFile", function(static, ...) {
  "[.](birdseed[.]data[.]txt)$";
}, static=TRUE)



setMethodS3("getReadArguments", "BroadBirdseedGenotypeTcgaDataFile", function(this, ..., colClassPatterns=c("*"="character", "Call$"="integer", "Confidence$"="double")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE)



setMethodS3("extractCalls", "BroadBirdseedGenotypeTcgaDataFile", function(this, ..., drop=TRUE) {
  colClassPatterns <- c("CompositeElement REF"="character", "Call$"="integer");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);
  idx <- match("CompositeElement REF", colnames(data));
  unitNames <- data[,idx];
  data <- data[,-idx,drop=FALSE];
  names <- names(data);
  pattern <- "(.*),(Call)$";
  sampleNames <- gsub(pattern, "\\1", names);
  nbrOfSamples <- length(sampleNames);

  # Coerce to a matrix  
  data <- as.matrix(data);
  rownames(data) <- unitNames;

  # A matrix? (probably never happens /HB 2009-08-23)
  if (drop && nbrOfSamples == 1) {
    data <- as.vector(data);
  }
  
  data;
})



setMethodS3("extractConfidenceScores", "BroadBirdseedGenotypeTcgaDataFile", function(this, ..., drop=TRUE) {
  colClassPatterns <- c("CompositeElement REF"="character", "Confidence$"="double");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);
  idx <- match("CompositeElement REF", colnames(data));
  unitNames <- data[,idx];
  data <- data[,-idx,drop=FALSE];
  names <- names(data);
  pattern <- "(.*),(Call)$";
  sampleNames <- gsub(pattern, "\\1", names);
  nbrOfSamples <- length(sampleNames);

  # Coerce to a matrix  
  data <- as.matrix(data);
  rownames(data) <- unitNames;

  # A matrix? (probably never happens /HB 2009-08-23)
  if (drop && nbrOfSamples == 1) {
    data <- as.vector(data);
  }
  
  data;
})




############################################################################
# HISTORY:
# 2009-10-25
# o Created.
############################################################################
