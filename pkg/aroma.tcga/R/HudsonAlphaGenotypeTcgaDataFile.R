setConstructorS3("HudsonAlphaGenotypeTcgaDataFile", function(...) {
  this <- extend(TcgaDataFile(...), "HudsonAlphaGenotypeTcgaDataFile");
  this;
})


setMethodS3("getExtensionPattern", "HudsonAlphaGenotypeTcgaDataFile", function(static, ...) {
  "[.](Genotypes[.]txt)$";
}, static=TRUE)



setMethodS3("getReadArguments", "HudsonAlphaGenotypeTcgaDataFile", function(this, ..., colClassPatterns=c("*"="character", "(Chr|Pos)"="integer", "genotype$"="character")) {
  NextMethod("getReadArguments", this, ..., colClassPatterns=colClassPatterns);
}, protected=TRUE)



setMethodS3("extractCalls", "HudsonAlphaGenotypeTcgaDataFile", function(this, sampleNames=NULL, ..., drop=TRUE) {
  # Argument 'sampleNames':
  if (!is.null(sampleNames)) {
    sampleNames <- Arguments$getCharacters(sampleNames);
  }


  colClassPatterns <- c("CompositeElement REF"="character");
  if (is.null(sampleNames)) {
    types <- c("genotype$"="character");
  } else {
    patterns <- sprintf("%s,genotype$", sampleNames);
    types <- rep("character", length(patterns));
    names(types) <- patterns;
  }
  colClassPatterns <- c(colClassPatterns, types);


  data <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);
  idx <- match("CompositeElement REF", colnames(data));
  unitNames <- data[,idx];
  data <- data[,-idx,drop=FALSE];
  names <- names(data);

  pattern <- "(.*),(genotype)$";
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
# 2009-12-05
# o Created.
############################################################################
