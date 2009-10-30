############################################################################
# Download missing files of existing TCGA data set directories
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");
require("R.menu") || throw("Package not loaded: R.menu");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scan for all TCGA archive directories
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "tcgaData";
rootPath <- Arguments$getReadablePath(rootPath);

paths <- TcgaDccData$findDataSets(rootPath=rootPath);
if (length(paths) == 0) {
  throw("No TCGA data set directories found.");
}
verbose && str(verbose, paths);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Check which ones don't have the file(s) of interest
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
filenames <- c("README.txt");
filenames <- NULL;
nbrOfFiles <- length(filenames);
if (nbrOfFiles > 0) {
  for (pp in seq(along=paths)) {
    path <- paths[pp];
    pathnames <- file.path(path, filenames);
    ok <- sapply(pathnames, FUN=isFile);
    if (all(ok)) {
      paths[pp] <- NA;
    }
  } # for (pp ...)
  verbose && str(verbose, paths);
  paths <- paths[!is.na(paths)];
}
verbose && str(verbose, paths);
if (length(paths) == 0) {
  throw("Nothing more to do.");
}

pattern <- sprintf("%s/", rootPath);
urlPaths <- gsub(pattern, "", paths);
verbose && str(verbose, urlPaths);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup all downloaders
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tddList <- lapply(urlPaths, FUN=TcgaDccDownloader);
username <- Sys.getenv("TCGAUSER");
pwd <- Sys.getenv("TCGAPWD");
tddList <- lapply(tddList, FUN=function(tdd) {
  setUsername(tdd, username);
  setPassword(tdd, pwd);
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tumorTypes <- sapply(tddList, FUN=getTumorType);
if (interactive()) {
  ans <- selectMenu(sort(unique(tumorTypes)));
  keep <- is.element(tumorTypes, ans);
  tddList <- tddList[keep];
} else {
  ans <- "gbm";
  keep <- is.element(tumorTypes, ans);
  tddList <- tddList[keep];
}
verbose && print(verbose, tddList);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select centers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
centers <- sapply(tddList, FUN=getCenter);
if (interactive()) {
  ans <- selectMenu(sort(unique(centers)));
  keep <- is.element(centers, ans);
  tddList <- tddList[keep];
} else {
  ans <- "broad.mit.edu";
  keep <- is.element(centers, ans);
  tddList <- tddList[keep];
}
verbose && print(verbose, tddList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop old data sets if newer exists
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- TcgaDccData$getDataSetPatterns()$dataset;
dataSets <- sapply(tddList, FUN=getDataSet);
dataSetNames <- gsub(pattern, "\\1_\\2.\\3", dataSets);
majorArchives <- as.integer(gsub(pattern, "\\5", dataSets));
labels <- sprintf("%s.%02d", dataSetNames, majorArchives);
o <- order(labels);
tddList <- tddList[o];

dataSets <- sapply(tddList, FUN=getDataSet);
uDataSets <- unique(dataSets);
res <- list();
for (kk in seq(along=dataSets)) {
  name <- dataSets[kk];
  res[[name]] <- c(res[[name]], tddList[kk]);
}
for (kk in seq(along=res)) {
  tddKK <- res[[kk]];
  if (length(tddKK) > 1) {
    versionsKK <- sapply(tddKK, FUN=function(tdd) {
      gsub(pattern, "\\6.\\7", getDataSet(tdd));
    });
    maxPos <- 1;
    for (jj in 2:length(versionsKK)) {
      if (compareVersion(versionsKK[jj], versionsKK[maxPos])) {
        maxPos <- jj;
      }
    }
    tddKK <- tddKK[[maxPos]];
  } else {
    tddKK <- tddKK[[1]];
  }
  res[[kk]] <- tddKK;
} # for (kk ...)
tddList <- res;
verbose && print(verbose, tddList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select data sets to download
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) {
  names <- names(tddList);
  names <- selectMenu(names, selected=rep(TRUE, length(names)));
  tddList <- tddList[names];
}
verbose && print(verbose, tddList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Download all core files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (pp in seq(along=tddList)) {
  tdd <- tddList[[pp]];
  urlPath <- getUrlPath(tdd);
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                        pp, urlPath, length(tddList)));

  verbose && enter(verbose, "Download core files");
  downloadCoreTextFile(tdd, verbose=less(verbose,10));
  verbose && exit(verbose);

  verbose && exit(verbose);
} # for (pp ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Download all missing files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tdd <- tddList[[1]];
patterns <- getKnownFilePatterns(tdd);
patterns <- unlist(patterns, use.names=FALSE);
if (interactive()) {
  patterns <- selectMenu(sort(patterns));
} else {
  patterns <- "[.](cel|CEL)$";
  patterns <- "[.]birdseed[.]data[.]txt$";
}

for (pp in seq(along=tddList)) {
  tdd <- tddList[[pp]];
	  urlPath <- getUrlPath(tdd);
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                        pp, urlPath, length(tddList)));

  verbose && enter(verbose, "Download files by patterns");
  for (kk in seq(along=patterns)) {
    pattern <- patterns[kk];
    verbose && enter(verbose, sprintf("Pattern #%d of %d", kk, length(patterns)));
    downloadFilesByPattern(tdd, pattern=pattern, order="random", 
                           onError="warning", verbose=less(verbose,10));
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  verbose && exit(verbose);
} # for (pp ...)


############################################################################
# HISTORY:
# 2009-10-23
# o Added interactive selection of download filename patterns.
# 2009-10-22
# o Now the MANIFEST.txt file is downloaded and parsed first of all.
# 2009-10-20
# o Created.
############################################################################
