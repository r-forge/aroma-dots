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
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPaths <- c("rawCnData", "cnData");
rootPath <- textMenu(rootPaths, value=TRUE);
rootPath <- Arguments$getReadablePath(rootPath);

# Find all data sets
paths <- list.files(path=rootPath, full.names=TRUE);
paths <- paths[sapply(paths, FUN=isDirectory)];

# Find TCGA data sets
tumorTypes <- c("OV", "GBM");
if (interactive()) {
  tumorType <- textMenu(tumorTypes, value=TRUE);
} else {
  tumorType <- tumorTypes[1];
}
tcgaPattern <- sprintf(".*/([a-z.]+)_(%s)[.]([^.]+)[.]([0-9.]+),.*,pairs", tumorType);
paths <- grep(tcgaPattern, paths, value=TRUE);

# Sanity check
if (length(paths) == 0) {
  throw("No TCGA data sets found.");
}

# Rules: For Broad, drop all non-5NN data sets
isSite <- (regexpr("broad.mit.edu", paths) != -1);
isType <- (regexpr("5NN", paths) != -1);
paths[isSite & !isType] <- NA;
paths <- paths[!is.na(paths)];

# Rules: For HudsonAlpha, drop all non-PCN data sets
isSite <- (regexpr("hudsonalpha.org", paths) != -1);
isType <- (regexpr("PCN", paths) != -1);
paths[isSite & !isType] <- NA;
paths <- paths[!is.na(paths)];

# Sanity check

if (length(paths) == 0) {
  throw("No TCGA data sets found after applying filter rules.");
}

verbose && cat(verbose, "Found ", length(paths), " TCGA data sets of interest.");
verbose && print(verbose, paths);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parse data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
# 2010-01-05
# o Added header information.
# 2009-10-05
# o Now it is possible to scan both rawCnData/ and cnData/.
# 2009-09-29
# o Created.
############################################################################
