############################################################################
# Author: Henrik Bengtsson, hb@stat.berkeley.edu
#
# Description:
# Generates one plot per sample and chromosome with smooth CN profiles
# for each of the sources before and after MSCN normalization.
# This script also do MSCN if not already done.
############################################################################
if (interactive()) savehistory();
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
writeCopyNumberProfiles <- function(dfList, setName, tags=NULL, chromosomes=1:24, clim=c(0,8), logBase=NULL, chr01Width=2000, height=400, ...) {
  # Argument 'height':
  height <- Arguments$getInteger(height);

  # Argument 'chromosomes':
  chromosomes <- Arguments$getIndices(chromosomes);

  df <- dfList[[1]];
  chipType <- getChipType(df);

  fullname <- paste(c(setName, tags), collapse=",");
  figPath <- filePath("reports", fullname, chipType);
  figPath <- Arguments$getWritablePath(figPath);

  chr01Length <- 247000001;
  for (cc in seq(along=chromosomes)) {
    chr <- chromosomes[cc];
    verbose && enter(verbose, sprintf("Chromosome #%d ('Chr%02d') of %d", 
                                           cc, chr, length(chromosomes)));
    figName <- sprintf("%s,chr%02d", sampleName, chr);
    filename <- sprintf("%s.png", figName);
    pathname <- filePath(figPath, filename);
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already plotted.");
      verbose && exit(verbose);
      next;
    }

    # Extract signals
    cnList <- lapply(dfList, FUN=extractRawCopyNumbers, chromosome=chr);
    verbose && print(verbose, cnList);

    # Merge and sort
    cn <- Reduce(append, cnList);
    cn <- sort(cn);
    verbose && print(verbose, cn);

    if (is.null(logBase)) {
      # Report on the non-log scale, i.e. C = 2*tumor/normal.
      cn$y <- 2*2^cn$y;
    }

    chrLength <- xMax(cn);
    width <- as.integer(chr01Width * chrLength / chr01Length);
    verbose && printf(verbose, "Image dimension: %dx%d pixels\n", height, width);

    devNew("png", pathname, width=width, height=height);
    par(mar=c(3,3,1,1)+0.1, mgp=c(1.8,0.5,0), cex=1.2);

    # Plot along chromosome
    plot(cn, col=cols[cn$id], pch=19, cex=0.7, xScale=1e-6, ylim=clim);
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, sprintf("Chr %02d", chr));
    labels <- gsub("/", "\n", sourceNames);
    legend("topleft", labels, col=cols, pch=19, cex=0.9, pt.cex=1.5, horiz=TRUE, bty="n");

    devDone();

    verbose && exit(verbose);
  } # for (chr ...)
} # writeCopyNumberProfiles()


timeStamp <- format(Sys.time(), format="%Y-%m-%d %k:%M:%S %Z");
pid <- sprintf("%s,%s", System$getHostname(), timeStamp);

analysisTag <- "HB20091001";

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
str(data);



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
filename <- sprintf("QC,smoothCNs,tracks,%s,jobs,queue.Rbin", analysisTag);
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
  colorMap <- c("mskcc.org/HG-CGH-1x1M"=1, "broad.mit.edu/GenomeWideSNP_6"=2, 
                "hms.harvard.edu/HG-CGH-244A"=3, "hudsonalpha.org/Human1M-Duo"=4,
                "hms.harvard.edu/HG-CGH-415K"=5, "mskcc.org/HG-CGH-244A"=6);
  allSourceNames <- names(colorMap);

  sourceNames <- dataKK$sourceName;
  verbose && print(verbose, sourceNames);

  idxs <- match(allSourceNames, sourceNames);
  idxs <- na.omit(idxs);
  targetDimension <- idxs[1];
  targetSourceName <- sourceNames[targetDimension];

  verbose && printf(verbose, "Target dimension: %d (%s)\n", 
                                 targetDimension, targetSourceName);
  verbose && print(verbose, dsList[[targetDimension]]);

  mscn <- MultiSourceCopyNumberNormalization(dsList, fitUgp=fitUgp, 
             targetDimension=targetDimension, align="byChromosome", 
                                         tags=c("*", analysisTag));
  verbose && print(verbose, mscn);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot smooth CNs before MSCN
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cols <- colorMap[sourceNames];
  verbose && print(verbose, cols);

  verbose && enter(verbose, "Plotting smooth CNs before MSCN");

  setName <- "TCGA,OV";

  dsSList <- getSmoothedDataSets(mscn, verbose=verbose);
  verbose && print(verbose, dsSList);

  dfList <- lapply(dsSList, getFile, 1);
  verbose && print(verbose, dfList);

  tags <- "raw";
#  writeCopyNumberProfiles(dfList, setName=setName, tags=tags, clim=c(0,8), logBase=NULL);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot smooth CNs after MSCN
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  verbose && enter(verbose, "Plotting smooth CNs fater MSCN (requires more smoothing)");

#  dsNList <- process(mscn, verbose=verbose);
#  verbose && print(verbose, dsNList);

  # Setup a fake MSCN (just for smoothing MSCN data)
#  mscnN <- MultiSourceCopyNumberNormalization(dsNList, fitUgp=fitUgp, 
#               targetDimension=targetDimension, align="byChromosome",
#                                           tags=c("*", analysisTag));
#  verbose && print(verbose, mscnN);

#  dsSList <- getSmoothedDataSets(mscnN, verbose=verbose);
#  verbose && print(verbose, dsSList);

#  dfList <- lapply(dsSList, getFile, 1);
#  verbose && print(verbose, dfList);

#  tags <- getTags(mscn, collapse=",");
#  writeCopyNumberProfiles(dfList, setName=setName, tags=tags, clim=c(0,8), logBase=NULL);

#  verbose && exit(verbose);

  verbose && exit(verbose);
} # for (kk ...)



############################################################################
# HISTORY:
# 2009-09-30
# o Created.
############################################################################
