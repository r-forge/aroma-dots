############################################################################
# Author: Henrik Bengtsson, hb@stat.berkeley.edu
#
# Description:
# Generates one plot per sample and chromosome with smooth CN profiles
# for each of the sources before and after MSCN normalization.
# This script also do MSCN if not already done.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
library("aroma.cn");
library("R.menu");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

analysisTag <- "HB20100124";

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

  fullname <- paste(setName, collapse=",");
  tags <- paste(tags, collapse=",");
  figPath <- filePath("reports", fullname, chipType, tags);
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

    ns <- sapply(cnList, nbrOfLoci);
    n <- ns[1];
    # Sanity check
    stopifnot(all(ns == n));

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
    plot(cn, col=cols[cn$id], pch=".", cex=3, ylim=clim);
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, sprintf("Chr %02d", chr));
    stext(side=4, pos=0, sprintf("n=%dx%d", length(dfList), n), cex=0.7);
    stext(side=4, pos=1, chipType, cex=0.7);
    labels <- gsub("/", "\n", sourceNames);
    legend("topleft", labels, col=cols, pch=19, cex=0.9, pt.cex=1.5, 
                      horiz=TRUE, bty="n");

    devDone();

    verbose && exit(verbose);
  } # for (chr ...)
} # writeCopyNumberProfiles()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all tumor-normal CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "smoothCnData";
rootPath <- Arguments$getReadablePath(rootPath, mustExist=TRUE);
patterns <- TcgaDccData$getDataSetPatterns();
pattern <- patterns$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);
verbose && cat(verbose, "Available data sets:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Raw or MSCN normalized?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
method <- "raw";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,.*", patterns$dataset);
types <- gsub(pattern, "\\2", dataSets);
uTypes <- sort(unique(types));

# User selected types
if (interactive()) {
  uTypes <- selectMenu(uTypes, selected=(length(uTypes) == 1), title="Select type:");
} else {
  if (!is.null(args$types)) {
    uTypes <- intersect(uTypes, args$types);
  }
}

keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);


keep <- is.element(types, uTypes);
dataSets <- dataSets[keep];
stopifnot(length(dataSets) > 0);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select by center
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- sprintf("^%s,.*", patterns$dataset);
centers <- gsub(pattern, "\\1", dataSets);
uCenters <- sort(unique(centers));
if (interactive()) {
  uCenters <- selectMenu(uCenters, selected=TRUE);
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
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  dataSets <- selectMenu(dataSets, selected=(length(dataSets) == 1)); 
} else {
  if (!is.null(args$dataSets)) {
    dataSets <- intersect(dataSets, args$dataSets);
  }
}


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
} else {
  if (!is.null(args$sampleNames)) {
    sampleNames <- intersect(sampleNames, args$sampleNames);
  }
}
stopifnot(length(sampleNames) > 0);

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


dataSets <- unname(sapply(dsList, FUN=getFullName));
verbose && cat(verbose, "Data sets to be processed:");
verbose && print(verbose, dataSets);
stopifnot(length(dataSets) > 0);

sampleNames <- lapply(dsList, FUN=getNames);
sampleNames <- sort(unique(unlist(sampleNames, use.names=FALSE)));
verbose && cat(verbose, "Samples to be processed:");
verbose && print(verbose, sampleNames);
stopifnot(length(sampleNames) > 0);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sanity checks
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ugpList <- lapply(dsList, getAromaUgpFile);
ugp <- ugpList[[1]];
for (kk in seq(along=ugpList)) {
  stopifnot(equals(ugpList[[kk]], ugp));
}
rm(ugpList);
verbose && print(verbose, ugp);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Precalculatings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chromosomes <- 1:24;
if (interactive()) {
  chromosomes <- selectMenu(chromosomes);
}


pattern <- sprintf("^%s,.*", patterns$dataset);
centers <- gsub(pattern, "\\1", dataSets);
chipTypes <- gsub(pattern, "\\3", dataSets);

pairs <- paste(centers, chipTypes, sep="/");

# OLD:
colorMap <- c(
  "mskcc.org/HG-CGH-1x1M"=1,
  "broad.mit.edu/GenomeWideSNP_6"=2, 
  "hms.harvard.edu/HG-CGH-244A"=3,
  "hudsonalpha.org/Human1M-Duo"=4,
  "hms.harvard.edu/HG-CGH-415K"=5,
  "mskcc.org/HG-CGH-244A"=6
);

colorMap <- c(
 "broad.mit.edu/Genome_Wide_SNP_6" = "red",
 "hms.harvard.edu/HG-CGH-244A" = "blue",
 "hudsonalpha.org/HumanHap550" = "orange",
 "mskcc.org/HG-CGH-244A" = "purple"
);

sourceColors <- colorMap[pairs];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process sample by sample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (ii in seq(along=sampleNames)) {
  sampleName <- sampleNames[ii];
  verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, sampleName, length(sampleNames)));

  dsListII <- lapply(dsList, FUN=function(ds) {
    extract(ds, sampleName, onMissing="NA");
  });

  pattern <- sprintf("^%s,.*", patterns$dataset);
  types <- gsub(pattern, "\\2", dataSets);
  tumorType <- types[1];
  # Sanity check
  stopifnot(all(types == tumorType));

  dfList <- lapply(dsList, FUN=function(ds) {
    dsT <- extract(ds, sampleName, onMissing="NA");
    getFile(dsT, 1);
  });

  keep <- sapply(dfList, isFile);
  dfList <- dfList[keep];
  verbose && print(verbose, dfList);

  # Sanity check
  stopifnot(length(dfList) > 0);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Annotiations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  sourceNames <- pairs[keep];
  cols <- sourceColors[keep];

  setName <- c("TCGA", tumorType, analysisTag);
  verbose && enter(verbose, "Creating PNGs");
  writeCopyNumberProfiles(dfList, setName=setName, tags=method, 
                                          chromosomes=chromosomes);
  verbose && exit(verbose);

  verbose && exit(verbose);
} # for (ii ...)



############################################################################
# HISTORY:
# 2009-09-30
# o Created.
############################################################################
