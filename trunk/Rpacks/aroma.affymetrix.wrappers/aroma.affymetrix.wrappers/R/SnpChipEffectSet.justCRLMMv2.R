setMethodS3("justCRLMMv2", "SnpChipEffectSet", function(this, units=NULL, ram=1, ..., balance=1.5, minLLRforCalls=c(5, 1, 5), recalibrate=FALSE, verbose=FALSE) {
  require("oligo") || throw("Package not loaded: oligo");

  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);
  chipType <- gsub(",monocell", "", chipType);

#  if (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1) {
#    throw("Cannot extract SnpQSet: Unsupported chip type: ", chipType);
#  }

  maleIndex <- c();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':  
  if (is.null(units)) {
    # Identify all SNP_A-* units (as what is returned by oligo)
    units <- indexOf(cdf, pattern="^SNP_A-");
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes by CRLMM");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(this);
  data <- data.frame(gender=rep("female", nbrOfArrays));
  phenoData <- new("AnnotatedDataFrame", data=data);

  pdPkgName <- oligo::cleanPlatformName(chipType);
  verbose && cat(verbose, "Platform Design (PD) package: ", pdPkgName);
 
  # Load target from PD package
  path <- system.file(package=pdPkgName);
  if (path == "") {
    throw("Cannot load HapMap reference target quantiles. Package not installed: ", pdPkgName);
  }

  verbose && enter(verbose, "Getting Platform Design Database");
  require(pdPkgName, character.only=TRUE) || throw("Package not loaded: ", pdPkgName);
  pdDB <- db(get(pdPkgName, mode="S4"));
  verbose && print(verbose, pdDB);
  verbose && exit(verbose);
 
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
  rm(pdPkgName);

  dim <- dim(crlmm$params$centers);
  hasQuartets <- (length(dim) == 3);

  verbose && enter(verbose, "Identifying all SNPs");
  # This can be inferred from the CDF as well.
  allSNPs <- dbGetQuery(pdDB, "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid")[[1]];
  verbose && str(verbose, allSNPs);

  # Sanity check
  if (length(allSNPs) != length(crlmm$hapmapCallIndex)) {
    throw("Internal error: The number of identified SNPs and the number of prior HapMap calls does not match: ", length(allSNPs), " != ", length(crlmm$hapmapCallIndex));
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying SNP on ChrX");
  # This can be inferred from the CDF and the UGP file.
  snpsOnChrX <- dbGetQuery(pdDB, "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' AND chrom = 'X'")[[1]];
  verbose && exit(verbose);
  rm(pdDB);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating parameter files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up output directory
  rootPath <- Arguments$getWritablePath("crlmmData");
  outPath <- filePath(rootPath, getFullName(this), chipType);
  outPath <- Arguments$getWritablePath(outPath);

  # Setting up output pathnames
  fullnames <- getFullNames(this);
  fullnames <- gsub(",chipEffects$", "", fullnames);
  filenames <- sprintf("%s,CRLMM.acu", fullnames);

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
      acu[,1] <- NA;
      verbose && exit(verbose);
    }
    acuList[[kk]] <- acu;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Processing units in chunk
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chunkSize <- ram * (500e3/nbrOfArrays);
  unitList <- splitInChunks(units, chunkSize=chunkSize);
  nbrOfChunks <- length(unitList);
  for (kk in seq(length=nbrOfChunks)) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", kk, nbrOfChunks));
    unitsKK <- unitList[[kk]];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsKK);

    verbose && enter(verbose, "Extracting data");
    if (hasQuartets) {
      eSet <- extractSnpQSet(this, units=unitsKK, sortUnits=FALSE, verbose=verbose);
    } else {
      eSet <- extractSnpCnvQSet(this, units=unitsKK, sortUnits=FALSE, verbose=verbose);
    }
    phenoData(eSet) <- phenoData;
    verbose && exit(verbose);

    unitNamesKK <- featureNames(eSet);
    verbose && cat(verbose, "Unit names:");
    verbose && str(verbose, unitNamesKK);

    verbose && enter(verbose, "Extract CRLMM priors");
    idxs <- match(unitNamesKK, allSNPs);
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

    indexX <- whichVector(is.element(unitNamesKK, snpsOnChrX));
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
    for (ii in seq(along=acuList)) {
      acu <- acuList[[ii]];
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, getName(acu), nbrOfArrays));
      # Calls
      acu[unitsKK,1] <- calls[,ii,drop=TRUE];
      # LLR
      acu[unitsKK,2] <- llr[,ii,drop=TRUE];
      # Distances
      distII <- dist[,ii,,,drop=TRUE];
      dim(distII) <- c(dim(distII)[1], prod(dim(distII)[2:3]));
      for (cc in 1:ncol(distII)) {
        acu[unitsKK,cc+2] <- distII[,cc];
      }
      rm(distII);
      verbose && exit(verbose);
    }
    verbose && exit(verbose);

    rm(unitsKK, calls, llr, dist);
    verbose && exit(verbose);
  } # for (kk ...)

  res <- AromaCrlmmBinarySet$fromFiles(outPath);

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2008-12-05
# o Created from justCRLMMv2() of oligo v1.7.3.
############################################################################
