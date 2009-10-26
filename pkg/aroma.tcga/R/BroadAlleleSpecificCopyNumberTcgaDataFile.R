setConstructorS3("BroadAlleleSpecificCopyNumberTcgaDataFile", function(...) {
  this <- extend(TcgaDataFile(...), "BroadAlleleSpecificCopyNumberTcgaDataFile");
  this;
})


setMethodS3("getExtensionPattern", "BroadAlleleSpecificCopyNumberTcgaDataFile", function(static, ...) {
  "[.](copynumber[.]byallele[.]data[.]txt)$";
}, static=TRUE)



setMethodS3("getReadArguments", "BroadAlleleSpecificCopyNumberTcgaDataFile", function(this, ..., colClassPatterns=c("*"="character", "(Chromosome|PhysicalPosition)"="integer", "(Signal_A|Signal_B)$"="double")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE)



setMethodS3("extractAlleleSpecificCopyNumbers", "BroadAlleleSpecificCopyNumberTcgaDataFile", function(this, ..., drop=TRUE) {
  colClassPatterns <- c("CompositeElement REF"="character", "(Signal_A|Signal_B)$"="double");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);
  idx <- match("CompositeElement REF", colnames(data));
  unitNames <- data[,idx];
  data <- data[,-idx];
  names <- names(data);
  # Sanity check (even number of columns)
  stopifnot(length(names) %% 2 == 0);

  pattern <- "(.*),(Signal_A|Signal_B)$";
  sampleNames <- gsub(pattern, "\\1", names);
  sampleNames <- sampleNames[seq(from=1, to=length(sampleNames), by=2)];
  nbrOfSamples <- length(sampleNames);

  alleles <- gsub(pattern, "\\2", names);
  alleles <- gsub("Signal_", "", alleles, fixed=TRUE);

  # Coerce to a matrix  
  data <- as.matrix(data);
  colnames(data) <- alleles;
  rownames(data) <- unitNames;

  # An array? (probably never happens /HB 2009-08-23)
  if (!drop || nbrOfSamples > 1) {
    dimnames <- list(unitNames, alleles[1:2], sampleNames);
    dim <- sapply(dimnames, FUN=length);
    dim(data) <- dim;
    dimnames(data) <- dimnames;
  }
  
  data;
})



setMethodS3("extractTotalAndFracB", "BroadAlleleSpecificCopyNumberTcgaDataFile", function(this, ..., drop=TRUE) {
  data <- extractAlleleSpecificCopyNumbers(this, ..., drop=FALSE);
  dim <- dim(data);
  dimnames <- dimnames(data);

  theta <- data[,1,,drop=FALSE] + data[,2,,drop=FALSE];
  beta <- data[,2,,drop=FALSE] / theta;
  rm(data);
  data <- c(theta, beta);
  rm(theta, beta);

  # Drop singleton dimensions?
  nbrOfSamples <- dim[3];
  if (drop && nbrOfSamples == 1) {
    dimnames <- dimnames[-3];
    dim <- dim[-3];
  }

  dim(data) <- dim;
  dimnames[[2]] <- c("total", "fracB");
  dimnames(data) <- dimnames;

  data;
})


############################################################################
# HISTORY:
# 2009-08-23
# o Created.
############################################################################
