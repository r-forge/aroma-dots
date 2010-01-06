############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/normalize/");
#
# Description:
# Applies MSCN normalization to a set of samples, which are specified by
# loaded data fram of pre-scanned samples.
############################################################################
if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

timeStamp <- format(Sys.time(), format="%Y-%m-%d %k:%M:%S %Z");
pid <- sprintf("%s,%s", System$getHostname(), timeStamp);

analysisTag <- "HB20091001";
analysisTag <- "HB20100105";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select set of samples to normalize
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- "^TCGA,[^,]+,CN,samples,([-0-9]+).Rbin$";
filenames <- list.files(pattern=pattern, path=".");
filenames <- sort(filenames, decreasing=TRUE);
if (interactive()) {
  filename <- textMenu(filenames, value=TRUE);
} else {
  filename <- filenames[1];
}
data <- loadObject(filename);
verbose && str(verbose, data);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify unique sample names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNames <- sort(unique(data[,"sampleName"]));
nbrOfSamples <- length(sampleNames);
verbose && cat(verbose, "Found ", nbrOfSamples, " unique samples:");
verbose && str(verbose, sampleNames);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get "job queue"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Grabbing samples from job queue");
filename <- sprintf("MSCN,%s,jobs,queue.Rbin", analysisTag);
pathname <- Arguments$getReadablePathname(filename, mustExist=FALSE);
if (isFile(pathname)) {
  # Load job queue
  jobs <- loadObject(pathname);
  chunks <- jobs$todo;
} else {
  # Split in 10 chunks;
  nbrOfChunks <- 10;
  chunkSize <- ceiling(nbrOfSamples/nbrOfChunks);
  chunks <- splitInChunks(seq(length=nbrOfSamples), chunkSize=chunkSize);
  jobs <- list(todo=chunks, done=NULL);
}

verbose && cat(verbose, "Jobs before:");
verbose && str(verbose, jobs);

# Work on chunk on top of queue
chunk <- chunks[[1]];
jobs$todo <- chunks[-1];
job <- list(chunk);
names(job) <- pid;
jobs$done <- c(jobs$done, job);

# Save job queue (for next host)
saveObject(jobs, file=pathname);

verbose && cat(verbose, "Jobs afterward:");
verbose && str(verbose, jobs);
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Samples to be processed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sampleNames <- sampleNames[chunk];
nbrOfSamples <- length(sampleNames);
verbose && cat(verbose, "Samples to be processed:");
verbose && print(verbose, sampleNames);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Target loci used for fitting the normalization curve
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitUgp <- AromaUgpFile$byChipType("GenericHuman", tags="200kb,HB20090503");
verbose && print(verbose, fitUgp);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process sample by sample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(length=nbrOfSamples)) {
  sampleName <- sampleNames[kk];
  verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", kk, sampleName, nbrOfSamples));

  keep <- (data[,"sampleName"] == sampleName);
  dataKK <- data[keep,,drop=FALSE];
  nbrOfSources <- nrow(dataKK);
  verbose && cat(verbose, "Number of sites that processed this sample: ", nbrOfSamples);
  if (nbrOfSources == 1) {
    verbose && cat(verbose, "Nothing to do. This sample was only process by one site: ", sampleName);
    verbose && exit(verbose);
    next;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Loading data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Loading data sets");
  rootPath <- "rawCnData";
  rootPath <- Arguments$getReadablePath(rootPath);
  pattern <- sprintf("^%s,.*,total[.]asb$", sampleName);
  dsList <- mapply(dataKK$dataSet, dataKK$chipType, FUN=function(dataSet, chipType) {
    path <- filePath(rootPath, dataSet, chipType);
    path <- Arguments$getReadablePath(path);
    ds <- AromaUnitTotalCnBinarySet$fromFiles(path, pattern=pattern);
    if (length(ds) != 1) {
      throw("Found more than one match data file: ", paste(getFullNames(ds), collapse=", "));
    }
    list(ds);
  });
  verbose && print(verbose, dsList);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Multi-source CN normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  idxs <- match(c("mskcc.org/HG-CGH-1x1M", "broad.mit.edu/GenomeWideSNP_6", 
                  "hms.harvard.edu/HG-CGH-244A", "hudsonalpha.org/Human1M-Duo",
                  "hms.harvard.edu/HG-CGH-415K", "mskcc.org/HG-CGH-244A"), 
                  dataKK$sourceName);
  idxs <- na.omit(idxs);
  targetDimension <- idxs[1];
  sourceNames <- sort(unique(dataKK$sourceName));
  verbose && printf(verbose, "Target dimension: %d (%s)\n", 
                    targetDimension, sourceNames[targetDimension]);
  verbose && print(verbose, dsList[[targetDimension]]);

  mscn <- MultiSourceCopyNumberNormalization(dsList, fitUgp=fitUgp, 
             targetDimension=targetDimension, align="byChromosome", 
                                         tags=c("*", analysisTag));
  verbose && print(verbose, mscn);
  
  dsNList <- process(mscn, verbose=verbose);
  verbose && print(verbose, dsNList);

  verbose && exit(verbose);
} # for (kk ...)



############################################################################
# HISTORY:
# 2009-09-30
# o Now 'targetDimension' is set according to a prioritized list.
# 2009-09-29
# o Updated to process sample by sample in a different way.
# 2009-05-21
# o Added verbose output.
# 2009-05-04
# o Updated.
# 2009-xx-xx
# o Created.
############################################################################
