############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/normalize/");
#
# Description:
# 1. Identifies all tumor-normal pairs.
# 2. Tabulates which data sources have each pair.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPaths <- c("rawCnData", "cnData");
if (interactive()) {
  rootPath <- textMenu(rootPaths, value=TRUE);
} else {
  rootPath <- rootPaths[1];
}
rootPath <- Arguments$getReadablePath(rootPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
patterns <- TcgaDccData$getDataSetPatterns();
pattern <- patterns$dataset;

dataSets <- list.files(path=rootPath, pattern=pattern);
verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s.*", patterns$dataset);
types <- gsub(pattern, "\\2", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1),
                       title="Select type:");
} else {
  if (!is.null(args$types)) {
    uTypes <- intersect(uTypes, args$types);
  }
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


# Human-style sorting
## try(dataSets <- gtools::mixedsort(dataSets), silent=TRUE);

tumorType <- uTypes[1];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select by center
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s.*", patterns$dataset);
centers <- gsub(pattern, "\\1", dataSets);
uCenters <- sort(unique(centers));
if (interactive()) {
  uCenters <- selectMenu(uCenters);
} else {
  if (!is.null(args$centers)) {
    uCenters <- intersect(uCenters, args$centers);
  }
}
verbose && cat(verbose, "TCGA center:");
verbose && print(verbose, uCenters);
dataSets <- dataSets[is.element(centers, uCenters)];

verbose && cat(verbose, "Filtered data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select preprocessing types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,", patterns$dataset);
types <- gsub(pattern, "", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1), title="Select preprocessing methods:");
} else {
  if (!is.null(args$methods)) {
    uTypes <- intersect(uTypes, args$methods);
  }
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);

# Human-style sorting
## try(dataSets <- gtools::mixedsort(dataSets), silent=TRUE);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSets <- selectMenu(dataSets, selected=(length(dataSets) == 1), title="Select dataset(s):");
}
stopifnot(length(dataSets) > 0);

verbose && cat(verbose, "Selected data sets:");
verbose && print(verbose, dataSets);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select samples to process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get all data sets
dsList <- lapply(dataSets, FUN=function(dataSet) {
  AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
});

sampleNames <- lapply(dsList, FUN=getNames);
sampleNames <- sort(unique(unlist(sampleNames, use.names=FALSE)));
stopifnot(length(sampleNames) > 0);

# User selection
if (interactive()) {
  sampleNames <- selectMenu(sampleNames, selected=(length(sampleNames) == 1));
}

# Keep only samples and data sets of interest
dsList <- lapply(dsList, extract, sampleNames, onMissing="drop");
ns <- sapply(dsList, length);
dsList <- dsList[ns > 0];
stopifnot(length(dsList) > 0);

dataSets <- unname(sapply(dsList, FUN=getFullName));
verbose && cat(verbose, "Data sets to be processed:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

verbose && cat(verbose, "Samples to be processed:");
verbose && print(verbose, sampleNames);
stopifnot(length(sampleNames) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parse data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
paths <- sapply(dsList, getPath);
paths <- sapply(paths, getParent);

verbose && enter(verbose, "Scanning for samples");
data <- NULL;
for (kk in seq(along=paths)) {
  path <- paths[kk];
  fullname <- basename(path);
  parts <- strsplit(fullname, split=",", fixed=TRUE)[[1]];
  name <- parts[1];
  tags <- paste(parts[-1], collapse=",");
  site <- gsub("^([a-z.]+).*", "\\1", name);
  chipType <- list.files(path=path);
  stopifnot(length(chipType) == 1);
  pathT <- filePath(path, chipType);
  filenames <- list.files(path=pathT);
  parts <- strsplit(filenames, split=",", fixed=TRUE);
  sampleNames <- sapply(parts, FUN=function(p) p[1]);
  dataKK <- cbind(site=site, chipType=chipType, sampleName=sampleNames,
                  sourceName=file.path(site, chipType), dataSet=fullname);
  data <- rbind(data, dataKK);
} # for (kk ...)
data <- as.data.frame(data, stringsAsFactors=FALSE);
verbose && str(verbose, data);
verbose && exit(verbose);

verbose && enter(verbose, "Saving data table of scanned samples");
filename <- sprintf("TCGA,%s,CN,samples,%s.Rbin", tumorType, Sys.Date());
pathname <- Arguments$getWritablePathname(filename, mustNotExist=FALSE);
saveObject(data, file=pathname);
verbose && cat(verbose, "Pathname: ", pathname);
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify unique sample names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNames <- sort(unique(data[,"sampleName"]));
nbrOfSamples <- length(sampleNames);
verbose && cat(verbose, "Found ", nbrOfSamples, " unique samples:");
verbose && str(verbose, sampleNames);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify unique source names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sourceNames <- sort(unique(data[,"sourceName"]));
nbrOfSources <- length(sourceNames);
verbose && cat(verbose, "Found ", nbrOfSources, " unique sources:");
verbose && print(verbose, sourceNames);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify which samples have been run on which sources
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
samples <- matrix(FALSE, nrow=nbrOfSamples, ncol=nbrOfSources);
rownames(samples) <- sampleNames;
colnames(samples) <- sourceNames;
for (name in colnames(samples)) {
  keys <- subset(data, sourceName == name, select=sampleName, drop=TRUE);
  samples[keys,name] <- TRUE;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Write summary to file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
t <- samples;
t[samples] <- "X";
t[!samples] <- "-";
colnames(t) <- base::LETTERS[1:ncol(t)];
t <- as.data.frame(t, stringsAsFactors=FALSE);
df <- data.frame(sampleName=rownames(samples));
df <- cbind(df, t);
df$counts <- rowSums(samples);
rownames(df) <- NULL;
names(sourceNames) <- colnames(t);
attr(df, "sourceNames") <- sourceNames;

sourceCounts <- colSums(samples);
names(sourceCounts) <- colnames(t);

# Reorder
o <- order(df$counts, decreasing=TRUE);
df <- df[o,,drop=FALSE];

dateStamp <- format(Sys.time(), format="%Y%m%d");
filename <- sprintf("TCGA,%s,%s,samples,%s.txt", tumorType, rootPath, dateStamp);
pathname <- Arguments$getWritablePathname(filename, mustNotExist=FALSE);

con <- file(filename, open="w");
cat(file=con, "#############################################################################\n");
cat(file=con, "# List of TCGA ", tumorType, " ", nbrOfSamples, 
              " unique samples with available copy number data from\n", sep="");
cat(file=con, "# each of the following sources:\n");
cat(file=con, "#\n");
labels <- sprintf("%s: %s [%d samples]\n", names(sourceNames), sourceNames, sourceCounts);
cat(file=con, paste("# ", labels, sep=""), sep="");
cat(file=con, "#\n");
cat(file=con, "# Do not edit this file; it was generated automatically.\n");
timeStamp <- format(Sys.time(), format="%Y-%m-%d %k:%M:%S %Z");
cat(file=con, "# Created on: ", timeStamp, "\n", sep="");
cat(file=con, "#############################################################################\n");
write.table(file=con, df, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t");
close(con);

verbose && cat(verbose, "Summary table written: ", pathname);

############################################################################
# HISTORY:
# 2010-01-20
# o Recreated using new selection methods.
# 2010-01-05
# o Added header information.
# 2009-10-05
# o Now it is possible to scan both rawCnData/ and cnData/.
# 2009-09-29
# o Created.
############################################################################
