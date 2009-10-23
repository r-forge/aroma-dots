setConstructorS3("BroadTotalCopyNumberTcgaDataFile", function(...) {
  this <- extend(TcgaDataFile(...), "BroadTotalCopyNumberTcgaDataFile");
  this;
})


setMethodS3("getFilenamePattern", "BroadTotalCopyNumberTcgaDataFile", function(static, ...) {
  ".*[.]after_5NN.copynumber[.]data[.]txt$";
}, static=TRUE)



setMethodS3("getReadArguments", "BroadTotalCopyNumberTcgaDataFile", function(this, ..., colClassPatterns=c("*"="character", "(Chromosome|PhysicalPosition)"="integer", "Signal$"="double")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE)



setMethodS3("extractTotalCopyNumbers", "BroadTotalCopyNumberTcgaDataFile", function(this, ..., drop=TRUE) {
  colClassPatterns <- c("CompositeElement REF"="character", "Signal$"="double");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);
  idx <- match("CompositeElement REF", colnames(data));
  unitNames <- data[,idx];
  data <- data[,-idx,drop=FALSE];
  names <- names(data);

  pattern <- "(.*),(Signal)$";
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



setMethodS3("extractTotalAndFracB", "BroadTotalCopyNumberTcgaDataFile", function(this, ..., drop=TRUE) {
  data <- extractTotalCopyNumbers(this, ..., drop=FALSE);
  dim <- dim(data);
  dimnames <- dimnames(data);

  # Setup (theta, beta)
  naValue <- as.double(NA);
  theta <- data;
  beta <- naValue;
  dim(beta) <- dim;
  data <- c(theta, beta);
  rm(theta, beta);

  # Restructure to (units, [theta,beta], arrays)
  dim <- c(dim, 2);
  dimnames[[3]] <- c("total", "fracB");
  dim(data) <- dim;
  dimnames(data) <- dimnames;
  data <- aperm(data, perm=c(1,3,2));

  # Drop singleton dimensions?
  nbrOfSamples <- dim[3];
  if (drop && nbrOfSamples == 1) {
    dimnames <- dimnames[-3];
    dim <- dim[-3];
    dim(data) <- dim;
    dimnames(data) <- dimnames;
  }

  data;
})


############################################################################
# HISTORY:
# 2009-09-24
# o Created from allele-specific ditto.
############################################################################
