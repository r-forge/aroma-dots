setConstructorS3("TotalCnSmoothing", function(dataSet=NULL, ..., targetUgp=NULL, sd=50e3, .reqSetClass="AromaTotalCnBinarySet") {
  if (!is.null(dataSet)) {
    # Argument 'targetUgp':
    if (!inherits(targetUgp, "AromaUgpFile")) {
      throw("Argument 'targetUgp' is not an AromaUgpFile: ", 
                                                    class(targetUgp)[1]);
    }
  }

  # Argument 'sd':
  sd <- Arguments$getDouble(sd, range=c(0,Inf));


  extend(AromaTransform(dataSet=dataSet, ..., .reqSetClass=.reqSetClass), "TotalCnSmoothing",
    .targetUgp = targetUgp,
    .sd = sd
  );
})


setMethodS3("getParameters", "TotalCnSmoothing", function(this, ...) {
  list(
    targetUgp = this$.targetUgp,
    sd = this$.sd
  );
}, private=TRUE);


setMethodS3("getAsteriskTags", "TotalCnSmoothing", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=NULL, ...);

  # Add class-specific tags

  params <- getParameters(this);
  # Parameter 'by'
  byTag <- grep("(b|kb|Mb)$", getTags(params$targetUgp), value=TRUE);
  if (length(byTag) > 0) {
    byTag <- sprintf("by=%s", byTag[1]);
  }

  # Parameter 'sd'
  sdTag <- sprintf("sd=%.1fkb", params$sd/1e3);

  tags <- c(tags, byTag, sdTag);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } 

  tags;
}, protected=TRUE) 



setMethodS3("getRootPath", "TotalCnSmoothing", function(this, ...) {
  "cnData";
}, private=TRUE)



setMethodS3("getTargetUgpFile", "TotalCnSmoothing", function(this, ...) {
  this$.targetUgp;
})

setMethodS3("getPath", "TotalCnSmoothing", function(this, create=TRUE, ...) {
  path <- NextMethod("getPath", this, create=FALSE, ...);
  path <- dirname(path);
  targetUgp <- getTargetUgpFile(this);
  chipType <- getChipType(targetUgp, fullname=FALSE);

  # The full path
  path <- filePath(path, chipType, expandLinks="any");

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this));
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  # Create path?
  if (create) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
}, private=TRUE)



setMethodS3("getTargetPositions", "TotalCnSmoothing", function(this, ..., force=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  params <- getParameters(this);
  targetUgp <- params$targetUgp;

  # For now, always use all chromosomes
  chromosomes <- NULL;

  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- getChromosomes(targetUgp);
  } else {
    chromosomes <- Arguments$getIndices(chromosomes);
  }
  nbrOfChromosomes <- length(chromosomes);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  res <- this$.targetPositions;
  if (!force && !is.null(res)) {
    return(res);
  }

  verbose && enter(verbose, "Identifying all target positions");
  verbose && cat(verbose, "Chromosomes:");
  verbose && print(verbose, chromosomes);

  verbose && print(verbose, targetUgp);

  res <- list();
  for (cc in chromosomes) {
    chrTag <- sprintf("Chr%02d", cc);
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", 
                                             cc, chrTag, nbrOfChromosomes));
    verbose && cat(verbose, "Target positions:");
    units <- getUnitsOnChromosome(targetUgp, chromosome=cc);
    xOut <- getPositions(targetUgp, units);
    verbose && str(verbose, xOut);
    res[[chrTag]] <- list(chromosome=cc, units=units, xOut=xOut);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  this$.targetPositions <- res;

  res;
}, protected=TRUE)




setMethodS3("process", "TotalCnSmoothing", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (isDone(this)) {
    dsOut <- getOutputDataSet(this);
    return(invisible(dsOut));
  }

  verbose && enter(verbose, "Smoothing copy-number towards set of target loci");

  params <- getParameters(this);

  verbose && print(verbose, "Input data set:");
  ds <- getInputDataSet(this);
  verbose && print(verbose, ds);

  verbose && enter(verbose, "Identifying all target positions");
  targetList <- getTargetPositions(this, ...);
  nbrOfChromosomes <- length(targetList);
  verbose && str(verbose, targetList);

  targetUgp <- params$targetUgp;
  platform <- getPlatform(targetUgp);
  chipType <- getChipType(targetUgp);
  nbrOfUnits <- nbrOfUnits(targetUgp);
  rm(targetUgp);
  verbose && cat(verbose, "Total number of target units:", nbrOfUnits);
  verbose && exit(verbose);

  nbrOfArrays <- length(ds);
  for (kk in seq(ds)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array %d ('%s') of %d", 
                                            kk, getName(df), nbrOfArrays));

    path <- getPath(this);
    filename <- getFilename(df, translate=TRUE);
    pathname <- Arguments$getReadablePathname(filename, path=path, 
                                                         mustExist=FALSE);
    className <- class(df)[1];
    clazz <- Class$forName(className);
    if (isFile(pathname)) {
      dfOut <- newInstance(clazz, filename=pathname);
      if (nbrOfUnits != nbrOfUnits(dfOut)) {
        throw("The number of units in existing output file does not match the number of units in the output file: ", nbrOfUnits, " != ", nbrOfUnits(dfOut));
      }
      verbose && cat(verbose, "Skipping already existing output file.");
      verbose && exit(verbose);
      next;
    }

    verbose && print(verbose, df);

    # Preallocate vector
    M <- rep(as.double(NA), nbrOfUnits);

    verbose && enter(verbose, "Reading and smoothing input data");
    for (cc in seq(along=targetList)) {
      target <- targetList[[cc]];
      chromosome <- target$chromosome;
      chrTag <- sprintf("Chr%02d", chromosome);
  
      verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", 
                                               cc, chrTag, nbrOfChromosomes));
      verbose && cat(verbose, "Extracting raw CNs:");
      rawCNs <- extractRawCopyNumbers(df, chromosome=chromosome, 
                                                  verbose=less(verbose, 10));
      verbose && print(verbose, rawCNs);
      verbose && summary(verbose, rawCNs);

      verbose && cat(verbose, "Smoothing CNs:");
      verbose && cat(verbose, "Target positions:");
      verbose && str(verbose, target$xOut);

      smoothCNs <- gaussianSmoothing(rawCNs, xOut=target$xOut, sd=params$sd, 
                                                  verbose=less(verbose, 20));
      verbose && print(verbose, smoothCNs);
      verbose && summary(verbose, smoothCNs);

      M[target$units] <- smoothCNs$cn;
      verbose && exit(verbose);
    } # for (cc ...)

    verbose && cat(verbose, "Smoothed CNs across all chromosomes:");
    verbose && str(verbose, M);
    verbose && summary(verbose, M);
    verbose && printf(verbose, "Missing values: %d (%.1f%%) out of %d\n", 
                   sum(is.na(M)), 100*sum(is.na(M))/nbrOfUnits, nbrOfUnits);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing smoothed data");

    verbose && enter(verbose, "Allocating ", className);
    verbose && cat(verbose, "Pathname: ", pathname);
    footer <- list(
      sourceDataFile=list(
        fullname=getFullName(df), 
        platform=getPlatform(df), 
        chipType=getChipType(df), 
        checksum=getChecksum(df)
      ), parameters=list(
        targetUgp=list(
          fullname=getFullName(params$targetUgp),
          platform=getPlatform(params$targetUgp),
          chipType=getChipType(params$targetUgp),
          checksum=getChecksum(params$targetUgp)
        ),
        sd=params$sd
      )
    );
    dfOut <- clazz$allocate(filename=pathname, nbrOfRows=nbrOfUnits, 
                            platform=platform, chipType=chipType, 
                            footer=footer, verbose=less(verbose, 50));
    verbose && exit(verbose);

    dfOut[,1] <- M;
    rm(M);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);

  dsOut <- getOutputDataSet(this);
  invisible(dsOut);
})



setMethodS3("getOutputFiles", "TotalCnSmoothing", function(this, ...) {
  NextMethod("getOutputFiles", pattern=".*[.]asb$", ...);
}, protected=TRUE) 




############################################################################
# HISTORY:
# 2008-05-23
# o Created.
############################################################################
