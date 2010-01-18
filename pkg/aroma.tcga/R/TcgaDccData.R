# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# TCGA data set and file name regular expressions
# [1] http://stackoverflow.com/questions/1313934/how-are-nested-capturing-groups-numbered-in-regular-expressions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setConstructorS3("TcgaDccData", function(...) {
  extend(Object(), "TcgaDccData");
})

setMethodS3("getDataSetPatterns", "TcgaDccData", function(static, version=c("3", "2"), ...) {
  # Argument 'version':
  version <- match.arg(version);

  patterns <- list();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # TCGA Data Set Patterns
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # TCGA center, e.g. broad.mit.edu, hudsonalpha.org
  patterns$center <- "[a-z.]+";
  # TCGA tumor tag, e.g. GBM, OV
  patterns$tumor <- "[A-Z]+";
  # TCGA platform, e.g. CGH-1x1M_G4447A, Genome_Wide_SNP_6
  patterns$platform <- "[A-Z][-_A-Za-z0-9]*";
  # TCGA archive version, e.g. 11.5.0, 1.0.0
  patterns$archive <- "([0-9]+)[.]([0-9]+)[.]([0-9]+)";
  patterns$level <- "(Level_[0-9])";
  
  patterns$datasetV2 <- with(patterns, {
    sprintf("(%s)_(%s)[.](%s)[.](%s)", 
             center, tumor, platform, archive);
  });

  patterns$datasetV3 <- with(patterns, {
    sprintf("(%s)_(%s)[.](%s)[.](%s)[.](%s)", 
             center, tumor, platform, level, archive);
  });

  if (version == 2) {
    patterns$dataset <-patterns$datasetV2;
  } else if (version == 3) {
    patterns$dataset <-patterns$datasetV3;
  }

  patterns;
})



setMethodS3("getDataSetPattern", "TcgaDccData", function(static, name, ...) {
  patterns <- getDataSetPatterns(static, ...);
  patterns[[name]];  
}, static=TRUE)


setMethodS3("findDataSets2", "TcgaDccData", function(static, pattern=NULL, rootPath="rawCnData", ...) {
  # Argument 'rootPath':
  rootPathF <- rootPath;
  rootPathF <- Arguments$getReadablePath(rootPathF);

  # Build pattern
  defaultPattern <- TcgaDccData$getDataSetPattern("dataset", ...);
  if (is.null(pattern)) {
    pattern <- defaultPattern;
  } else {
    pattern <- sprintf(pattern, defaultPattern);
    pattern <- Arguments$getRegularExpression(pattern);
  }

  # Scan for directories matching the regular expression for TCGA DCC data sets
  paths <- list.files(pattern=pattern, path=rootPathF, full.names=TRUE);
  if (length(paths) > 0) {
    paths <- paths[sapply(paths, FUN=isDirectory)];
  }
  dataSets <- basename(paths);

  # For each data set, identify the chip type
  chipTypes <- character(length(paths));
  for (kk in seq(along=paths)) {
    path <- paths[kk];
    pathsKK <- list.files(path=path, full.names=TRUE);
    pathsKK <- pathsKK[sapply(pathsKK, FUN=isDirectory)];
    stopifnot(length(pathsKK) == 1);
    pathKK <- pathsKK[1];
    chipType <- basename(pathKK);
    chipTypes[kk] <- chipType;
    paths[kk] <- pathKK;
  } # for (kk ...)

  rootPaths <- rep(rootPath, times=length(paths));
  paths <- file.path(rootPath, dataSets, chipTypes);
  data <- data.frame(rootPath=rootPaths, dataSet=dataSets, 
                    chipType=chipTypes, path=paths, stringsAsFactors=FALSE);

  data;
}, static=TRUE) # findDataSets2()



############################################################################
# HISTORY:
# 2010-01-17
# o Added data set patterns for new DCC v3 archive names (with Level_N).
# 2009-10-20
# o Added new static findDataSets().
# o "Old" findDataSets() was renamed to findDataSets2().
# 2009-10-03
# o Added static findDataSets().
# 2009-10-02
# o Created.
############################################################################
