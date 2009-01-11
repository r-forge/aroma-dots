setMethodS3("exportAromaUnitSignalBinaryFileList", "SnpChipEffectFile", function(this, whats=c("total", "fracB"), fullname=getFullName(this), dataSet=NULL, path=NULL, ..., overwrite=FALSE, drop=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'fullname':
  fullname <- Arguments$getCharacter(fullname);

  # Argument 'what':
  whats <- match.arg(whats, several.ok=TRUE);

  # Arguments 'dataSet':
  if (is.null(dataSet)) {
    dataSet <- basename(getParent(getPath(this)));
  }
  dataSet <- Arguments$getCharacter(dataSet);

  # Arguments 'path':
  if (is.null(path)) {
    rootPath <- "cnData";
    chipType <- getChipType(this, fullname=FALSE);
    path <- filePath(rootPath, dataSet, chipType);
  }
  path <- Arguments$getWritablePath(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && cat(verbose, "Whats: ", paste(whats, collapse=", "));

  cdf <- getCdf(this); 
  nbrOfUnits <- nbrOfUnits(cdf);
  platform <- getPlatform(this);
  chipType <- getChipType(this, fullname=FALSE);

  footer <- list(
    srcFile=list(
      srcDataSet=dataSet,
      srcChipType=getChipType(this),
      srcFullName=getFullName(this),
      srcChecksum=getChecksum(this)
    )
  );

  data <- NULL;

  asbList <- list();
  for (what in whats) {
    # Identify output class
    if (what == "total") {
      signalClass <- AromaUnitTotalCnBinaryFile;
    } else if (is.element(what, c("fracB", "freqB"))) {
      signalClass <- AromaUnitFracBCnBinaryFile;
    }
  
    verbose && enter(verbose, "Exporting ", class(this)[1], " as an ", getName(signalClass));
    verbose && cat(verbose, "Signal: ", what);
  
    # Generate output filename
    filename <- sprintf("%s,%s.asb", fullname, what);
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE); 
    verbose && cat(verbose, "Output pathname: ", pathname);
    if (isFile(pathname)) {
      if (!overwrite) {
        verbose && cat(verbose, "Output file already exists. Return that instead.");
        asb <- signalClass$fromFile(pathname);
        asbList[[what]] <- asb;
        verbose && exit(verbose);
        next;
      }
    }
  
    verbose && cat(verbose, "File footer:");
    verbose && str(verbose, footer);
  
    # Reading data
    if (is.null(data)) {
      verbose && enter(verbose, "Reading data");
      data <- extractTotalAndFracB(this, verbose=less(verbose, 5));
      verbose && str(verbose, data);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Allocating output file");
    asb <- signalClass$allocate(filename=pathname, path=NULL, nbrOfRows=nbrOfUnits, platform=platform, chipType=chipType, footer=footer, overwrite=overwrite, verbose=less(verbose, 25));
    verbose && exit(verbose);
  
    verbose && enter(verbose, "Writing data");
    asb[,1] <- data[,what, drop=TRUE];
    verbose && exit(verbose);
    
    verbose && exit(verbose);

    asbList[[what]] <- asb;
  } # for (what ...)
  names(asbList) <- whats;
  rm(data);

  if (drop && length(asbList) == 1) {
    asbList <- asbList[[1]];
  }

  invisible(asbList);
}, protected=TRUE)



setMethodS3("exportAromaUnitSignalBinarySetList", "SnpChipEffectSet", function(this, whats=c("total", "fracB"), ..., drop=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'whats':
  whats <- match.arg(whats, several.ok=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  signalClassList <- lapply(whats, FUN=function(what) {
    if (what == "total") {
      signalClass <- AromaUnitTotalCnBinarySet;
    } else if (is.element(what, c("fracB", "freqB"))) {
      signalClass <- AromaUnitFracBCnBinarySet;
    }
    signalClass;
  });

  names <- paste(sapply(signalClassList, FUN=getName), collapse=" and ");
  verbose && enter(verbose, "Exporting ", class(this)[1], " as ", names);

  dataSetName <- getFullName(this);
  chipType <- NULL;
  for (kk in seq(this)) {
    cf <- getFile(this, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(cf), length(this)));
    asbList <- exportAromaUnitSignalBinaryFileList(cf, whats=whats, dataSet=dataSetName, ..., drop=FALSE, verbose=less(verbose, 1));
    if (is.null(chipType)) {
      chipType <- getChipType(asbList[[1]]);
    }
    rm(asbList);
    verbose && exit(verbose);
  }
  verbose && exit(verbose);
  

  assList <- lapply(signalClassList, function(signalClass) {
    verbose && enter(verbose, "Setting up the ", getName(signalClass));
    ass <- NULL;
    tryCatch({
      ass <- signalClass$byName(dataSetName, chipType=chipType);
      verbose && print(verbose, ass);
    }, error = function(ex) {
    })
    verbose && exit(verbose);
    ass;
  });
  names(assList) <- whats;

  assList <- assList[!sapply(assList, is.null)];

  if (drop && length(assList) == 1) {
    assList <- assList[[1]];
  }

  invisible(assList);
}, protected=TRUE)


setMethodS3("getAromaTotalCnBinarySet", "SnpChipEffectSet", function(this, ...) {
  exportAromaUnitSignalBinarySetList(this, whats="total", ...);
})

setMethodS3("getAromaFracBCnBinarySet", "SnpChipEffectSet", function(this, ...) {
  exportAromaUnitSignalBinarySetList(this, whats="fracB", ...);
})


setMethodS3("getAromaUnitSignalBinarySetList", "SnpChipEffectSet", function(this, whats=c("total", "fracB"), ...) {
  exportAromaUnitSignalBinarySetList(this, whats=whats, ...);
})


setMethodS3("getTotalAndFracBSets", "SnpChipEffectSet", function(this, whats=c("total", "fracB"), ...) {
  exportAromaUnitSignalBinarySetList(this, whats=whats, ...);
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Backward compatibility
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setMethodS3("getAromaFreqBCnBinarySet", "default", function(this, ...) {
  getAromaFracBCnBinarySet(this, ...);
}, protected=TRUE)

setMethodS3("getTotalAndFreqBSets", "default", function(this, ...) {
  getTotalAndFracBSets(this, ...);
}, protected=TRUE)

setMethodS3("exportAromaSignalBinaryFileList", "default", function(this, ...) {
  exportAromaUnitSignalBinaryFileList(this, ...);
}, protected=TRUE, deprecated=TRUE);

setMethodS3("exportAromaUnitSignalBinarySetList", "default", function(this, ...) {
  exportAromaUnitSignalBinarySetList(this, ...);
}, protected=TRUE, deprecated=TRUE);


############################################################################
# HISTORY:
# 2008-09-10
# o Updated to be compatible with new aroma.core.
# 2008-07-30
# o Added getTotalAndFreqBSets() which is a more convenient name.
# 2008-06-25
# o Created.
############################################################################ 
