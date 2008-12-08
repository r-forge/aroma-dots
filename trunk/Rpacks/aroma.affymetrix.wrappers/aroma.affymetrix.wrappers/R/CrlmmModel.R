setConstructorS3("CrlmmModel", function(dataSet=NULL, balance=1.5, minLLRforCalls=c(5, 1, 5), recalibrate=FALSE, flavor="v2", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    className <- "SnpChipEffectSet";
    if (!inherits(dataSet, className))
      throw("Argument 'dataSet' is not an ", className, ": ",
                                                       class(dataSet)[1]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Sanity check
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    chipType <- getChipType(dataSet);
    chipType <- gsub(",monocell", "", chipType);

    # For now, only allow know SNP chip types. /HB 2008-12-07
    if (regexpr("^Mapping(10|50|250)K_.*$", chipType) != -1) {
    } else if (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1) {
    } else {
      throw("Cannot fit CRLMM model: Unsupported/unsafe chip type: ", chipType);
    }
  }

  # Argument 'balance':
  balance <- Arguments$getDouble(balance, range=c(0.00001, 1e6));

  # Argument 'minLLRforCalls':
  minLLRforCalls <- Arguments$getDoubles(minLLRforCalls, length=3, 
                                                         range=c(0, 1e6));

  # Argument 'recalibrate':
  recalibrate <- Arguments$getLogical(recalibrate);

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  extend(Model(dataSet=dataSet, ...), "CrlmmModel",
    balance = balance, 
    minLLRforCalls = minLLRforCalls,
    recalibrate = recalibrate,
    flavor = flavor
  )
})



setMethodS3("getParameterSet", "CrlmmModel", function(this, ...) {
  params <- NextMethod("getParameterSet", this, ...);
  params$balance <- this$balance;
  params$minLLRforCalls <- this$minLLRforCalls;
  params$recalibrate <- this$recalibrate;
  params$flavor <- this$flavor;
  params;
}, private=TRUE) 




setMethodS3("getChipType", "CrlmmModel", function(this, ...) {
  ds <- getDataSet(this);
  chipType <- getChipType(cdf);
  chipType <- gsub(",monocell", "", chipType);
  chipType;
})



setMethodS3("getCallSet", "CrlmmModel", function(this, ..., verbose=FALSE) {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  ces <- getDataSet(this);
  cdf <- getCdf(ces);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating parameter files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up output directory
  rootPath <- Arguments$getWritablePath("crlmmData");
  outPath <- filePath(rootPath, getFullName(ces), chipType);
  outPath <- Arguments$getWritablePath(outPath);

  # Setting up output pathnames
  fullnames <- getFullNames(ces);
  fullnames <- gsub(",chipEffects$", "", fullnames);
  filenames <- sprintf("%s,CRLMM.acu", fullnames);

  nbrOfArrays <- nbrOfArrays(ces);
  nbrOfUnits <- nbrOfUnits(cdf);
  platform <- getPlatform(cdf);

  verbose && enter(verbose, "Retrieving genotype data set");
  acuList <- list();
  for (kk in seq(along=filenames)) {
    filename <- filenames[kk];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, filename, nbrOfArrays));
    pathname <- filePath(outPath, filename);
    if (isFile(pathname)) {
      acu <- AromaCrlmmBinaryFile(pathname);
    } else {
      verbose && enter(verbose, "Allocating new file");
      acu <- AromaCrlmmBinaryFile$allocate(filename=pathname, platform=platform, chipType=chipType, nbrOfRows=nbrOfUnits, verbose=log);
      verbose && exit(verbose);
    }
    acuList[[kk]] <- acu;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  res <- AromaCrlmmBinarySet$fromFiles(outPath);

  res;
})


setMethodS3("getUnitsToFit", "CrlmmModel", function(this, safe=TRUE, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (safe) {
    # Identify which units CRLMM have fitted
    crlmmSNPs <- getCrlmmSNPs(this, verbose=less(verbose, 1));
    units <- crlmmSNPs;
  } else {
    ds <- getDataSet(this);
    cdf <- getCdf(ds);
    # Identify al SNP_A-* units.
    units <- indexOf(cdf, pattern="^SNP_A-");

##    # Identify all genotyping units (SNPs)
##    unitTypes <- getUnitTypes(cdf, verbose=verbose);
##    snpUnits <- whichVector(unitTypes == 2);
  }

  units;
})


setMethodS3("findUnitsTodo", "CrlmmModel", function(this, units=NULL, safe=TRUE, ...) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  unitsToFit <- getUnitsToFit(this, ...);
  verbose && str(verbose, unitsToFit);

  if (safe) {
    unitsTodo <- unitsToFit;
  } else {
    # Find all units with missing-value (NA) calls 
    callSet <- getCallSet(this);
    unitsWithNAs <- findUnitsTodo(callSet, ...);
    unitsTodo <- intersect(unitsToFit, unitsWithNAs);
    verbose && str(verbose, unitsTodo);
  }

  if (!is.null(units)) {
    unitsTodo <- intersect(units, unitsTodo);
    verbose && str(verbose, unitsTodo);
  }

  unitsTodo;
})


setMethodS3("fit", "CrlmmModel", function(this, units="remaining", force=FALSE, ram=1, ..., verbose=FALSE) {
  require("oligo") || throw("Package not loaded: oligo");

  maleIndex <- c();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':  
  doRemaining <- FALSE;
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  } else if (identical(units, "remaining")) {
    doRemaining <- TRUE;
  } else {
    throw("Unknown mode of argument 'units': ", mode(units)); 
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes by CRLMM");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get algorithm parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this);
  balance <- params$balance;
  minLLRforCalls <- params$minLLRforCalls;
  recalibrate <- params$recalibrate;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying units to process
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units to process");
  if (is.null(units)) {
    unitsToFit <- getUnitsToFit(this, verbose=less(verbose,1));
    units <- unitsToFit;
    rm(unitsToFit);
  } else if (doRemaining) {
    verbose && enter(verbose, "Identifying non-estimated units")
    units <- findUnitsTodo(this, safe=FALSE, verbose=less(verbose));
    verbose && str(verbose, units);
    verbose && exit(verbose); 
  } else {
    unitsToFit <- getUnitsToFit(this, verbose=less(verbose,1));
    units <- intersect(units, unitsToFit);
    rm(unitsToFit);
    # Fit only unique units
    units <- unique(units); 
  }
  nbrOfUnits <- length(units); 
  verbose && str(verbose, units);
  verbose && printf(verbose, "Getting model fit for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  if (!doRemaining) {
    if (force) {
      verbose && printf(verbose, "All of these are forced to be fitted.\n");
    } else {
      units <- findUnitsTodo(this, units=units, safe=FALSE, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits);
    }
  } 
  verbose && exit(verbose);

  # Nothing to do?
  if (nbrOfUnits == 0) {
    verbose && exit(verbose);
    return(invisible(NULL));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setup");
  crlmm <- getCrlmmPriors(this, verbose=less(verbose,1));

  ces <- getDataSet(this);
  nbrOfArrays <- nbrOfArrays(ces);
  data <- data.frame(gender=rep("female", nbrOfArrays));
  phenoData <- new("AnnotatedDataFrame", data=data);


  allSNPs <- getCrlmmSNPs(this, verbose=less(verbose, 1));
  verbose && str(verbose, allSNPs);
  snpsOnChrX <- getCrlmmSNPsOnChrX(this, verbose=less(verbose,1));
  verbose && str(verbose, snpsOnChrX);

  # Sanity check
  if (length(allSNPs) != length(crlmm$hapmapCallIndex)) {
    throw("Internal error: The number of identified SNPs and the number of prior HapMap calls does not match: ", length(allSNPs), " != ", length(crlmm$hapmapCallIndex));
  }

  dim <- dim(crlmm$params$centers);
  hasQuartets <- (length(dim) == 3);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get result set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callSet <- getCallSet(this, verbose=less(verbose,1));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Processing units in chunk
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chunkSize <- ram * (500e3/nbrOfArrays);
  unitList <- splitInChunks(units, chunkSize=chunkSize);
  nbrOfChunks <- length(unitList);
  count <- 1;
  while (length(unitList) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks));
    unitsChunk <- unitList[[1]];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsChunk);
    unitList <- unitList[-1];

    verbose && enter(verbose, "Extracting data");
    if (hasQuartets) {
      eSet <- extractSnpQSet(ces, units=unitsChunk, sortUnits=FALSE, verbose=verbose);
    } else {
      eSet <- extractSnpCnvQSet(ces, units=unitsChunk, sortUnits=FALSE, verbose=verbose);
    }
    phenoData(eSet) <- phenoData;
    verbose && exit(verbose);

    unitNamesKK <- featureNames(eSet);
    verbose && cat(verbose, "Unit names:");
    verbose && str(verbose, unitNamesKK);

    verbose && enter(verbose, "Extract CRLMM priors");
    idxs <- match(unitNamesKK, names(allSNPs));
    hapmapCallIndex <- crlmm$hapmapCallIndex[idxs];
    params <- crlmm$params;
    if (hasQuartets) {
      params$centers <- params$centers[idxs,,];
      params$scales <- params$scales[idxs,,];
    } else {
      params$centers <- params$centers[idxs,];
      params$scales <- params$scales[idxs,];
    }
    params$N <- params$N[idxs,];
    rm(idxs);
    verbose && cat(verbose, "CRLMM prior parameters (estimated from HapMap):");
    verbose && str(verbose, params);
    verbose && exit(verbose);

    verbose && enter(verbose, "Fitting SNP mixtures");
    correction <- oligo:::fitAffySnpMixture(eSet, verbose=as.logical(verbose));
    verbose && str(verbose, correction);
    verbose && exit(verbose);

    verbose && enter(verbose, "Genotype calling");
    naValue <- as.integer(NA);
    calls <- matrix(naValue, nrow=nrow(eSet), ncol=ncol(eSet));
    index <- whichVector(!hapmapCallIndex);
    if (length(index) > 0) {
      verbose && enter(verbose, "Initial SNP calling");
      verbose && str(verbose, index);
      calls[index,] <- oligo:::getInitialAffySnpCalls(correction, index, sqsClass=class(eSet), verbose=as.logical(verbose));
      verbose && exit(verbose);
    }
    verbose && exit(verbose);

    verbose && enter(verbose, "Estimate genotype regions");
    rparams <- oligo:::getAffySnpGenotypeRegionParams(eSet, calls, correction$fs, subset=index, sqsClass=class(eSet), verbose=as.logical(verbose));
    verbose && exit(verbose);

    priors <- crlmm$priors;
    if (hasQuartets) {
      verbose && enter(verbose, "Updating sense & antisense genotype regions");
      eSet1 <- eSet[,1];
      M <- getM(eSet1);
      dimnames(M) <- NULL;
      M <- M[,1,,drop=TRUE];
      oneStrand <- integer(nrow(M));
      for (ss in 1:2) {
        oneStrand[is.na(M[,ss])] <- ss;
      }
      rparams <- oligo:::updateAffySnpParams(rparams, priors, oneStrand, verbose=as.logical(verbose));
      params  <- oligo:::replaceAffySnpParams(params, rparams, index);
      dist <- oligo:::getAffySnpDistance(eSet, params, correction$fs);
      dist[,,-2,] <- balance*dist[,,-2,];
      verbose && exit(verbose);
    } else {
      verbose && enter(verbose, "Updating genotype regions");
      rparams <- oligo:::updateAffySnpParams(rparams, priors, verbose=as.logical(verbose));
      params  <- oligo:::replaceAffySnpParams(params, rparams, index);
      dist <- oligo:::getAffySnpDistance(eSet, params, correction$fs);
      dist[,,-2] <- balance*dist[,,-2];
      verbose && exit(verbose);
    }
    rm(params, index);

    indexX <- whichVector(is.element(unitNamesKK, names(snpsOnChrX)));
    calls <- oligo:::getAffySnpCalls(dist, indexX, maleIndex, sqsClass=class(eSet), verbose=as.logical(verbose));
    llr <- oligo:::getAffySnpConfidence(dist, calls, indexX, maleIndex, verbose=as.logical(verbose));
    
    if (recalibrate) {
      verbose && enter(verbose, "Recalibrating");
      naValue <- as.integer(NA);
      for (gg in 1:3) {
        calls[calls == gg & llr < minLLRforCalls[gg]] <- naValue;
      }
      rm(llr);
      calls[,(correction$snr < 3.675)] <- naValue;

      rparams <- oligo:::getAffySnpGenotypeRegionParams(eSet, calls, correction$fs, verbose=as.logical(verbose));
      rm(calls);

      rparams <- oligo:::updateAffySnpParams(rparams, priors, oneStrand);
      dist <- oligo:::getAffySnpDistance(eSet, rparams, correction$fs, verbose=as.logical(verbose));
      if (hasQuartets) {
        dist[,,-2,] <- balance*dist[,,-2,];
      } else {
        dist[,,-2] <- balance*dist[,,-2];
      }
      rm(oneStrand);
      calls <- oligo:::getAffySnpCalls(dist, indexX, maleIndex, verbose=as.logical(verbose));
      llr <- oligo:::getAffySnpConfidence(dist, calls, indexX, maleIndex, verbose=as.logical(verbose));
      verbose && exit(verbose);
    } # if (recalibrate)

    # Clean up
    rm(eSet, indexX,  correction, rparams, priors);

    verbose && enter(verbose, "Estimated genotype parameters");
    verbose && cat(verbose, "dist:");
    verbose && str(verbose, dist);
    verbose && cat(verbose, "calls:");
    verbose && str(verbose, calls);
    verbose && cat(verbose, "llr:");
    verbose && str(verbose, llr);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing genotype parameters");
    for (ii in seq(callSet)) {
      acu <- getFile(callSet, ii);
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, getName(acu), nbrOfArrays));
      # Calls
      acu[unitsChunk,1] <- calls[,ii,drop=TRUE];
      # LLR
      acu[unitsChunk,2] <- llr[,ii,drop=TRUE];
      # Distances
      distII <- dist[,ii,,,drop=TRUE];
      dim(distII) <- c(dim(distII)[1], prod(dim(distII)[2:3]));
      for (cc in 1:ncol(distII)) {
        acu[unitsChunk,cc+2] <- distII[,cc];
      }
      rm(distII, acu);
      verbose && exit(verbose);
    }
    verbose && exit(verbose);

    # Next chunk
    count <- count + 1;
    rm(unitsChunk, calls, llr, dist);
    verbose && exit(verbose);
  } # while(length(unitsList) > 0)
  rm(callSet);

  verbose && exit(verbose);

  # Return fitted units
  invisible(units);
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Plaform Design package dependent code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getCrlmmPriors", "CrlmmModel", function(this, ..., verbose=FALSE) {
  chipType <- getChipType(this);
  pdPkgName <- oligo::cleanPlatformName(chipType);
  verbose && cat(verbose, "Platform Design (PD) package: ", pdPkgName);
 
  # Load target from PD package
  path <- system.file(package=pdPkgName);
  if (path == "") {
    throw("Cannot load HapMap reference target quantiles. Package not installed: ", pdPkgName);
  }

  verbose && enter(verbose, "Loading CRLMM priors etc");
  path <- file.path(path, "extdata");
  path <- Arguments$getReadablePath(path);
  filename <- sprintf("%sCrlmmInfo.rda", pdPkgName);
  pathname <- Arguments$getReadablePathname(filename, path=path);
  verbose && cat(verbose, "Pathname: ", pathname);
  key <- sprintf("%sCrlmm", pdPkgName);
  crlmm <- loadToEnv(pathname)[[key]];
  verbose && capture(verbose, ll(envir=crlmm));
  verbose && exit(verbose);

  crlmm;
}, protected=TRUE);


setMethodS3("getPlatformDesignDB", "CrlmmModel", function(this, ..., verbose=FALSE) {
  verbose && enter(verbose, "Getting Platform Design Database");
  chipType <- getChipType(this);
  verbose && cat(verbose, "Chip type: ", chipType);
  pdPkgName <- oligo::cleanPlatformName(chipType);
  verbose && cat(verbose, "Plaform Design package: ", pdPkgName);
  require(pdPkgName, character.only=TRUE) || throw("Package not loaded: ", pdPkgName);
  pdDB <- db(get(pdPkgName, mode="S4"));
  verbose && print(verbose, pdDB);
  verbose && exit(verbose);
  pdDB;
}, private=TRUE)


setMethodS3("getCrlmmSNPs", "CrlmmModel", function(this, ..., verbose=FALSE) {
  # This can be inferred from the CDF as well.

  verbose && enter(verbose, "Identifying SNP on ChrX");
  pdDB <- getPlatformDesignDB(this, verbose=less(verbose,1));
  res <- dbGetQuery(pdDB, "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid")[[1]];
  verbose && str(verbose, res);

  ds <- getDataSet(this);
  cdf <- getCdf(this);
  units <- indexOf(cdf, names=res);
  names(units) <- res;

  verbose && str(verbose, units);

  verbose && exit(verbose);

  units;
}, private=TRUE)


setMethodS3("getCrlmmSNPsOnChrX", "CrlmmModel", function(this, ..., verbose=FALSE) {
  # This can be inferred from the CDF and the UGP file.

  verbose && enter(verbose, "Identifying all SNPs");
  pdDB <- getPlatformDesignDB(this, verbose=less(verbose,1));
  res <- dbGetQuery(pdDB, "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' AND chrom = 'X'")[[1]];
  verbose && str(verbose, res);

  ds <- getDataSet(this);
  cdf <- getCdf(this);
  units <- indexOf(cdf, names=res);
  names(units) <- res;

  verbose && str(verbose, units);

  verbose && exit(verbose);

  units;
}, private=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Plaform Design package dependent code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

############################################################################
# HISTORY:
# 2008-12-08
# o Now setup is much more like fit() for ProbeLevelModel.
# 2008-12-07
# o Starting to make use of AromaCrlmmBinarySet.
# o Created CrlmmModel from justCRLMMv2().
# 2008-12-05
# o Created from justCRLMMv2() of oligo v1.7.3.
############################################################################
