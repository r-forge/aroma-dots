dataSetName <- "HapMap270,6.0,CEU,founders"
chipType <- "GenomeWideSNP_6"; chipTypeTag="Full";

dataPath <- Arguments$getReadablePath("rocData", dataSetName);

if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType, tags=chipTypeTag);
  print(cdf);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get annotation data for units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("units23")) {
  pathname <- filePath(dataPath, sprintf("%s,units23.RData", getChipType(cdf)));
  units23 <- loadObject(pathname);
  str(units23);
}

if (!exists("units23excl")) {
  # Names of the units not in GTC
  filename <- "GenomeWideSNP_6,SNPsCNsMissingInGTC.txt";
  pathname <- Arguments$getReadablePathname(filename, path=dataPath);
  unExcl <- readLines(pathname);
  units23excl <- indexOf(cdf, names=unExcl);
  rm(unExcl);
  stopifnot(!anyMissing(units23excl));
  str(units23excl);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get genomic positions of the units on ChrX
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("positions", mode="numeric")) {
  gi <- getGenomeInformation(cdf);
  print(gi);
  positions <- getPositions(gi, units=units23);
  str(positions);
  rm(gi);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Identify known CNV regions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("cnvData")) {
  # Names of the units not in GTC
  cnvPath <- filePath("cnvData/RedonR_etal_2006,CNV_Project/NCBI36");
  filename <- "Redon_CNVs.gff";
  pathname <- Arguments$getReadablePathname(filename, path=cnvPath);
  columnNames <- c("chr", "study", "index", "start", "end", 
                                        "V6", "strand", "V8", "V9");
  cnv <- TabularTextFile(pathname, skip=2, columnNames=columnNames);
  colClassPatterns <- c("chr"="character", "(start|end)"="integer");
  data <- readDataFrame(cnv, colClassPatterns=colClassPatterns);
  data$chr <- gsub("^chr", "", data$chr);
  data$chr <- gsub("X", "23", data$chr);
  data$chr <- gsub("Y", "24", data$chr);
  data$chr <- as.integer(data$chr);
  cnvData <- data;
  rm(data);
  str(cnvData);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Exclude some regions on ChrX
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("regionsExcludedOnChr23")) {
  # Exclude known CNV regions with a safety margin
  regions <- subset(cnvData, chr == 23)[,c("start", "end")];
  regions <- as.matrix(regions);
  margin <- 100e3;
  regions[,"start"] <- regions[,"start"]-margin;
  regions[,"end"] <- regions[,"end"]+margin;

  # Exclude pseudo copy-neutral regions on ChrX
  par1 <- as.integer(c(0,2692881));
  par2 <- as.integer(c(154494747,154824264));
  regions <- rbind(regions, par1);
  regions <- rbind(regions, par2);

  # Sort regions
  regions <- regions[order(regions[,"start"]),];
  rownames(regions) <- NULL;

  # Merge regions (assuming already ordered)
  regions2 <- matrix(as.integer(0), nrow=0, ncol=2);
  colnames(regions2) <- colnames(regions);
  currRegion <- regions[1,,drop=FALSE];
  for (kk in seq(from=2, to=nrow(regions))) {
    nextRegion <- regions[kk,];

    # Does the next region overlap with the current one?
    if (nextRegion[1] <= currRegion[2]) {
      # Does it stretch beyond the current one?
      if (nextRegion[2] > currRegion[2]) {
        # Then merge the two
        currRegion[2] <- nextRegion[2];
        nextRegion <- NULL;
      } else {
        # Drop the next region because it is fully 
        # included in the current one.
        nextRegion <- NULL;
      }
    } else {
      # The next and current regions are disjoint.
      regions2 <- rbind(regions2, currRegion);
      currRegion <- nextRegion;
    }
  } # for (kk ...)
  regions2 <- rbind(regions2, currRegion);
  rownames(regions2) <- NULL;

  regionsExcludedOnChr23 <- regions2;
  rm(regions, regions2, nextRegion, currRegion);
}

# Keep only units with positions outside the exclude regions
exclA <- rep(FALSE, length=length(positions));
for (kk in 1:nrow(regionsExcludedOnChr23)) {
  region <- regionsExcludedOnChr23[kk,];
  isInside <- (region[1] < positions & positions < region[2]);
  exclA[isInside] <- TRUE;
}
cat("Units excluded because they are within a CNV+-100kb:\n");
print(summary(exclA));
##    Mode   FALSE    TRUE
## logical   71192   16012   


# Exclude units not processed by Affymetrix' GTC
exclB <- (units23 %in% units23excl);
cat("Units excluded because they are processed by Affymetrix' GTC:\n");
print(summary(exclB));
##    Mode   FALSE    TRUE
## logical   85712    1492

cat("Total set of units excluded:\n");
excl <- (exclA | exclB);
print(summary(excl));
##    Mode   FALSE    TRUE
## logical   70107   17097

cat("Units kept for this study:\n");
keep <- !excl;
print(summary(keep));
keep <- which(keep);
str(keep);

# Sanity check
if (length(keep)/length(positions) < 1/4) {
  throw(sprintf("Internal error. Too few loci: %.2f%%", 100*sum(keep)));
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Filter (positions, units23) accordingly
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
positions <- positions[keep];
units23 <- units23[keep];
unitsToKeep <- keep;
rm(keep, excl, exclA, exclB);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Reorder (positions, units23, unitsToKeep) by positions
# (this will allow up to smooth data directly by rows)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
reorder <- TRUE;
if (reorder && is.unsorted(positions)) {
  o <- order(positions);
  positions <- positions[o];
  units23 <- units23[o];
  unitsToKeep <- unitsToKeep[o];
  rm(o);
}
nbrOfUnits <- length(units23);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get unit types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("unitTypes", mode="integer")) {
  # A bit slow, right now (so infer from unit names instead)
  ## unitTypes <- getUnitTypes(cdf, units=units23);

  # Default is CN probes
  unitTypes <- rep(as.integer(5), length(units23));
  unitNames <- getUnitNames(cdf, units=units23);
  unitTypes[grep("SNP", unitNames)] <- as.integer(2);
  rm(unitNames);
  attr(unitTypes, "map") <- c(snp=as.integer(2), cn=as.integer(5));
  str(unitTypes);
}
knownUnitTypes <- c("all", "snp", "cn");


#####################################################################
# NOTE: Only (positions, unitTypes) & unitsToKeep is used later.
#####################################################################




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get annotation data for samples
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("sampleNames", mode="character")) {
  pathname <- filePath(dataPath, "sampleNames.RData");
  sampleNames <- loadObject(pathname);
  print(sampleNames);
  rm(n23);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Exclude bad samples
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
samplesToExclude <- "NA12145";
sampleNames <- sampleNames[!sampleNames %in% samplesToExclude];
rm(samplesToExclude);
nbrOfSamples <- length(sampleNames);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get known ChrX ploidy, i.e. gender
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("n23")) {
  pathname <- "annotationData/samples/HapMap270.saf";
  saf <- SampleAnnotationFile$fromFile(pathname);
  print(saf);
#  df <- readDataFrame(saf, colClassPattern=c("name"="character", "gender"="character"));  # TO IMPLEMENT
  df <- readDataFrame(saf);

  # Keep only samples to be studied (in same order)
  df <- df[match(sampleNames, df[,"name"]),];

  # Infer the ChrX ploidy
  n23 <- as.integer(c(male=1, female=2)[df[,"gender"]]);
  names(n23) <- df[,"name"];
  print(n23);
  rm(df, saf);
  gc <- gc();
}



############################################################################
# HISTORY: 
# 2008-07-18
# o Created from old CRMA and CRMA6 scripts.
############################################################################
