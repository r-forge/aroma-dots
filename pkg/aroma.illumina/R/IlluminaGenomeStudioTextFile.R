setConstructorS3("IlluminaGenomeStudioTextFile", function(...) {
  extend(TabularTextFile(...), "IlluminaGenomeStudioTextFile");
})


setMethodS3("getSampleNames", "IlluminaGenomeStudioTextFile", function(this, ...) {
  names <- getColumnNames(this, ...);
  pattern <- "[.]GType$";
  names <- grep(pattern, names, value=TRUE);
  names <- gsub(pattern, "", names);
  names;
})

setMethodS3("nbrOfSamples", "IlluminaGenomeStudioTextFile", function(this, ...) {
  length(getSampleNames(this, ...));
})



setMethodS3("readUnitNames", "IlluminaGenomeStudioTextFile", function(this, ...) {
  patterns <- c(Index="integer", Name="character");
  data <- readDataFrame(this, colClassPatterns=patterns, ...);
  unitNames <- data$Name;
  unitNames;
}) # readUnitNames()


setMethodS3("readUnitGenomePositions", "IlluminaGenomeStudioTextFile", function(this, ...) {
  patterns <- c(Index="integer", Name="character", Chr="character", Position="integer");
  data <- readDataFrame(this, colClassPatterns=patterns, ...);
  unitNames <- data$Name;
  chrs <- data$Chr;
  chrs[chrs == "X"] <- 23L;
  chrs[chrs == "Y"] <- 24L;
  chrs[chrs == "M"] <- 25L;
  gp <- data.frame(chromosome=chrs, position=data$Position, stringsAsFactors=FALSE);
  rownames(gp) <- unitNames;
  gp;
}) # readUnitGenomePositions()


setMethodS3("readSampleData", "IlluminaGenomeStudioTextFile", function(this, sample, patterns=c("GType"="character", "Score"="double", "Log R Ratio"="double", "B Allele Freq"="double"), ...) {
  # Argument 'sample':
  sampleNames <- getSampleNames(this);
  if (is.numeric(sample)) {
    sample <- Arguments$getIndex(sample, max=length(sampleNames));
    sample <- sampleNames[sample];
  } else if (!is.element(sample, sampleNames)) {
    throw("Unknown sample name: ", sample);    
  }


  # Read sample columns
  fields <- names(patterns);
  names(patterns) <- sprintf("^%s.%s", sample, names(patterns));
  patterns <- c(Name="character", patterns);

  data <- readDataFrame(this, colClassPatterns=patterns, ...);
  unitNames <- data$Name;

  data <- data[-1L];
  rownames(data) <- unitNames;
  pattern <- sprintf("^%s.", sample);
  colnames(data) <- gsub(pattern, "", colnames(data));

  data;
})


############################################################################
# HISTORY:
# 2012-02-13
# o Created.
############################################################################
