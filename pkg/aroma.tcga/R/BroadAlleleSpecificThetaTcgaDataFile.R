setConstructorS3("BroadAlleleSpecificThetaTcgaDataFile", function(...) {
  this <- extend(TcgaDataFile(...), "BroadAlleleSpecificThetaTcgaDataFile");
  this;
})


setMethodS3("getFilenamePattern", "BroadAlleleSpecificThetaTcgaDataFile", function(static, ...) {
  "[.]ismpolish[.]data[.]txt$";
}, static=TRUE)



setMethodS3("getReadArguments", "BroadAlleleSpecificThetaTcgaDataFile", function(this, ..., colClassPatterns=c("*"="character", "(Chromosome|PhysicalPosition)"="integer", "(Signal_A|Signal_B)$"="double")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE)



setMethodS3("extractTheta", "BroadAlleleSpecificThetaTcgaDataFile", function(this, ..., drop=TRUE) {
  colClassPatterns <- c("CompositeElement REF"="character", "Signal$"="double");
  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);

  # Sanity check
  stopifnot(ncol(data) == 2);

  sampleNames <- gsub(",Signal$", "", colnames(data)[-1]);

  unitNames <- data[,1,drop=TRUE];
  values <- data[,2,drop=TRUE];
  rm(data);

  # Identify alleles (if SNPs)
  alleles <- rep(as.character(NA), length(unitNames));
  len <- nchar(unitNames);
  tail <- substr(unitNames, start=len-1, stop=len);
  for (aa in c("A", "B")) {
    idxs <- whichVector(tail == sprintf("-%s", aa));
    alleles[idxs] <- aa;
    unitNames[idxs] <- substr(unitNames[idxs], start=1, stop=len[idxs]-2);
    rm(idxs);
  }

  # Identify the unique vector of unit names
  uniqueUnitNames <- unique(unitNames);
  nbrOfUnits <- length(uniqueUnitNames);

  # Map unit names to this vector
  units <- match(unitNames, uniqueUnitNames);
  rm(unitNames);

  # Allocate and populate Jx2 theta matrix
  naValue <- as.double(NA);
  data <- matrix(naValue, nrow=nbrOfUnits, ncol=2);
  idxs <- whichVector(is.na(alleles));
  data[units[idxs],1] <- values[idxs];
  idxs <- whichVector(alleles == "A");
  data[units[idxs],1] <- values[idxs];
  idxs <- whichVector(alleles == "B");
  data[units[idxs],2] <- values[idxs];
  rownames(data) <- uniqueUnitNames;
  rm(units, uniqueUnitNames, idxs);

  # Return an array?
  if (!drop) {
    dimnames <- dimnames(data);
    dimnames[[3]] <- sampleNames;
    dim <- c(dim(data), length(sampleNames));
    dim(data) <- dim;
    dimnames(data) <- dimnames;
  }

  data;
})



setMethodS3("extractTotalAndFracB", "BroadAlleleSpecificThetaTcgaDataFile", function(this, ..., drop=TRUE) {
  data <- extractTheta(this, ..., drop=FALSE);

  # Add data[,2,] to the total data[,1,], if non-NA.
  ok <- whichVector(!is.na(data[,2,]));
  data[ok,1,] <- data[ok,1,] + data[ok,2,];
  rm(ok);

  # Calculate allele B fraction
  data[,2,] <- data[,2,] / data[,1,];

  # Drop singleton dimensions?
  nbrOfSamples <- dim[3];
  if (drop && nbrOfSamples == 1) {
    dimnames <- dimnames[-3];
    dim <- dim[-3];
  }

  dimnames(data)[[2]] <- c("total", "fracB");

  data;
})


############################################################################
# HISTORY:
# 2009-08-26
# o Implemented extractTheta() and extractTotalAndFracB() correctly.
# 2009-08-23
# o Created.
############################################################################
