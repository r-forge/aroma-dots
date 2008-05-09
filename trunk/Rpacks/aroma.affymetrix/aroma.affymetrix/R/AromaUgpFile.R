setConstructorS3("AromaUgpFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUgpFile",
    .chromosomes=NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})

setMethodS3("getChromosomes", "AromaUgpFile", function(this, force=FALSE, .chromosomes=NULL, ...) {
  chromosomes <- this$.chromosomes;
  if (force || is.null(chromosomes)) {
    chromosomes <- .chromosomes;
    if (is.null(chromosomes))
      chromosomes <- this[,1,drop=TRUE];
    chromosomes <- unique(chromosomes);
    chromosomes <- chromosomes[!is.na(chromosomes)];
    chromosomes <- sort(chromosomes);
    this$.chromosomes <- chromosomes;
  }
  chromosomes;
})

setMethodS3("getFilenameExtension", "AromaUgpFile", function(static, ...) {
  "ugp";
}, static=TRUE)

setMethodS3("getColumnNames", "AromaUgpFile", function(this, ...) {
  c("chromosome", "position");
})

setMethodS3("readDataFrame", "AromaUgpFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  data <- NextMethod("readDataFrame", this, ..., verbose=less(verbose));

  verbose && enter(verbose, "Converting zeros to NAs");
  # Interpret zeros as NAs
  if (ncol(data) > 0) {
    for (cc in seq(length=ncol(data))) {
      nas <- (!is.na(data[,cc]) & (data[,cc] == 0));
      data[nas,cc] <- NA;
    }
  }
  verbose && exit(verbose);

  data;
})

setMethodS3("getGenomeVersion", "AromaUgpFile", function(this, ...) {
  tags <- getTags(this, ...);
  tags <- grep("^hg", tags, value=TRUE);
  tags;
}, protected=TRUE)


setMethodS3("allocateFromCdf", "AromaUgpFile", function(static, ...) {
  # NextMethod() not supported here.
  allocateFromCdf.AromaUnitTabularBinaryFile(static, ..., types=rep("integer",2), sizes=c(1,4));
}, static=TRUE)



setMethodS3("getUnitsAt", "AromaUgpFile", function(this, chromosomes, region=NULL, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument '...':
  args <- list(...);
  if ("chromosome" %in% names(args)) {
    chromosomes <- args[["chromosome"]];
  }


  # Stratify by chromosome
  data <- this[,1,drop=TRUE];

  # Update known chromosomes, if not already done.  
  allChromosomes <- getChromosomes(this, .chromosomes=data);

  keep <- !is.na(data) & (data %in% chromosomes);
  idxs <- which(keep);

  if (!is.null(region)) {
    data <- this[idxs,2,drop=TRUE];
    keep <- !is.na(data);
    keep <- keep & (region[1] <= data & data <= region[2]);
    idxs <- idxs[keep];
  }
  
  idxs;
}, protected=TRUE)




setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUgpFile", function(this, csv, shift="auto", onReplicates=c("median", "mean", "overwrite"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  if (!inherits(csv, "AffymetrixNetAffxCsvFile")) {
    throw("Argument 'csv' is not an AffymetrixNetAffxCsvFile: ", class(csv)[1]);
  }

  # Argument 'shift':
  if (!identical(shift, "auto"))
    shift <- Arguments$getInteger(shift);

  # Argument 'onReplicates':
  onReplicates <- match.arg(onReplicates);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Importing (unit name, chromosome, position) data from ", class(csv)[1]);

  # Query CDF
  cdf <- getCdf(this);
  cdfUnitNames <- getUnitNames(cdf);

  # Read data
  data <- readDataUnitChromosomePosition(csv, ..., verbose=less(verbose));
  importNames <- attr(data, "importNames");

  # Map to CDF unit names
  cdfUnits <- match(data[,1], cdfUnitNames);

  # Exclude units that are not in the CDF
  keep <- which(!is.na(cdfUnits));
  cdfUnits <- cdfUnits[keep];
  if (length(cdfUnits) == 0) {
    warning("None of the imported unit names match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }
  data <- data[keep,2:3,drop=FALSE];
  importNames <- importNames[2:3];

  # Garbage collect
  rm(keep);
  gc <- gc();


  # Shift positions?
  if (identical(shift, "auto")) {
    shift <- 0;
    if (any(regexpr("[sS]tart", importNames) != -1))
      shift <- 13;
  }
  if (shift != 0) {
    verbose && printf(verbose, "Shifting positions %d steps.", shift);
    data[,2] <- data[,2] + as.integer(shift);
  }
 
  # Replicated positions per unit?
  dups <- which(duplicated(cdfUnits));
  if (length(dups) > 0) {
    verbose && enter(verbose, "Detected units with replicated positions");
    dupUnits <- unique(cdfUnits[dups]);
    nDupUnits <- length(dupUnits);
    verbose && cat(verbose, "Number of units with replicated positions: ", 
                                                               nDupUnits);

    if (onReplicates %in% c("median", "mean")) {
      verbose && enter(verbose, "Calculate average positions for those (assuming they are on the same chromosome)");
      if (onReplicates == "mean") {
        avgFcn <- median;
      } else {
        avgFcn <- mean;
      }
      for (kk in seq(along=dupUnits)) {
        if (kk %% 500 == 0)
         verbose && printf(verbose, "%d, ", kk);
        dupUnit <- dupUnits[kk];
        # Identify position
        units <- which(cdfUnits == dupUnit);
        # Average position
        avgPos <- median(data[units,2], na.rm=TRUE);
        avgPos <- round(avgPos);
        # Update (can we update just units[1]?)
        data[units,2] <- avgPos;
      }
      verbose && cat(verbose, kk);
      verbose && exit(verbose);
      verbose && enter(verbose, "Remove the extraneous cases");
      data <- data[-dups,,drop=FALSE];
      cdfUnits <- cdfUnits[-dups];
      verbose && exit(verbose);
      rm(dupUnits, units);
      verbose && str(verbose, cdfUnits);
      warning("The positions for ", nDupUnits, " units were calculated as the average of replicated positions, since that was what was available on file.");
    } else {
      verbose && cat(verbose, "Ignored replicated positions. The last one written will be the one available.");
    }
    verbose && exit(verbose);
  }
  rm(dups);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  this[cdfUnits,1] <- data[,1];
  this[cdfUnits,2] <- data[,2];

  rm(data);
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  verbose && exit(verbose);

  invisible(cdfUnits);
})



setMethodS3("importFromAffymetrixTabularFile", "AromaUgpFile", function(this, src, colClassPatterns=c("*"="NULL", "^probeSetID$"="character", "^chromosome$"="character", "^(physicalPosition|position)$"="character"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!inherits(src, "AffymetrixTabularFile")) {
    throw("Argument 'src' is not an AffymetrixTabularFile: ", class(src)[1]);
  }

  units <- importFromGenericTabularFile(this, src=src, 
            colClassPatterns=colClassPatterns, camelCaseNames=TRUE, ...);

  invisible(units);
}, protected=TRUE);



setMethodS3("importFromGenericTabularFile", "AromaUgpFile", function(this, src, colClassPatterns=c("*"="NULL", "^Probe Set ID$"="character", "^Chromosome$"="character", "^Physical Position$"="character"), colOrder=NULL, shift=0, con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!inherits(src, "GenericTabularFile")) {
    throw("Argument 'src' is not an GenericTabularFile: ", class(src)[1]);
  }

  # Argument 'colOrder':
  if (!is.null(colOrder)) {
    colOrder <- Arguments$getIndices(colOrder, length=3);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Importing (unitName, chromosome, position) from ", class(src)[1], " file");

  data <- readDataFrame(src, colClassPatterns=colClassPatterns, ..., verbose=less(verbose));

  # Rearrange columns (optional)
  if (!is.null(colOrder))
    data <- data[,colOrder,drop=FALSE];

  # Map to CDF unit names
  cdf <- getCdf(this);
  cdfUnitNames <- getUnitNames(cdf);
  cdfUnits <- match(data[,1], cdfUnitNames);

  # Exclude units that are not in the CDF
  keep <- which(!is.na(cdfUnits));
  cdfUnits <- cdfUnits[keep];
  if (length(cdfUnits) == 0) {
    warning("None of the imported unit names match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(src), "'), being imported?");
  }

  # Assume 'chromosome' is in 2nd column, and 'position' in 3rd.
  data <- data[keep,2:3,drop=FALSE];

  # Garbage collection
  rm(keep);
  gc <- gc();

  # Convert chromosome strings to integers
  if (!is.integer(data[,1])) {
    map <- c(X=23, Y=24, Z=25);
    for (kk in seq(along=map)) {
      data[,1] <- gsub(names(map)[kk], map[kk], data[,1]);
    }
    suppressWarnings({
      data[,1] <- as.integer(data[,1]);
    })
    gc <- gc();
  }

  # Convert positions to integers
  if (!is.integer(data[,2])) {
    suppressWarnings({
      data[,2] <- as.integer(data[,2]);
    })
    gc <- gc();
  }

  # Shift positions?
  if (shift != 0) {
    verbose && printf(verbose, "Shifting positions %d steps.", shift);
    data[,2] <- data[,2] + as.integer(shift);
  }

  # Update to file
  this[cdfUnits,1] <- data[,1];
  this[cdfUnits,2] <- data[,2];

  verbose && exit(verbose);

  invisible(cdfUnits);
}, protected=TRUE);



setMethodS3("importFromGenomeInformation", "AromaUgpFile", function(this, gi, ..., verbose=FALSE) {
  if (!inherits(gi, "GenomeInformation")) {
    throw("Argument 'gi' is not a GenomeInformation object: ", class(gi)[1]);
  }

  # AD HOC patch, since units==NULL does not work./HB 2007-03-03
  units <- seq_len(nbrOfUnits(gi));
  data <- getData(gi, units=units, fields=c("chromosome", "physicalPosition"));

  chr <- data[,"chromosome"];
  if (is.character(chr)) {
    chr[chr == "X"] <- 23;
    chr[chr == "Y"] <- 24;
    chr[chr == "Z"] <- 25;
    suppressWarnings({
      chr <- as.integer(chr);
    })
  }
  
  pos <- data[,"physicalPosition"];
  suppressWarnings({
    pos <- as.integer(pos);
  })

  this[,1] <- chr;
  this[,2] <- pos;

  # A best guess of what was imported
  units <- units[!(is.na(chr) & is.na(pos))];

  invisible(units);
})


############################################################################
# HISTORY:
# 2008-04-29
# o BUG FIX: Name clash in getUnitsAt() after new argument 'chromosomes'.
# 2008-04-17
# o Renamed argument 'chromosome' of getUnitsAt() of AromaUgpFile to 
#   'chromosomes'.  This was done in order to make it consistent with 
#   getUnitsOnChromosome() of GenomeInformation. Thanks Tim Keighley at 
#   CSIRO for pointing this out.
# 2008-04-14
# o Renamed readData() to readDataFrame() for AromaTabularBinaryFile.
# 2007-09-16
# o Now importFromAffymetrixNetAffxCsvFile() averages positions if multiple
#   positions were available for a particular unit.
# o Now importFromGenomeInformation() tries to return units imported and
#   not all.  This is still a best guess, but still more informative than
#   before.
# 2007-09-14
# o Added getChromosomes(), which caches results in memory.
# o Added importFromAffymetrixTabularFile() to AromaUgpFile.
# 2007-09-13
# o Removed createFromGenomeInformation().
# o Updated AromaUgpFile according to changes in super class.
# 2007-09-10
# o Added importFromAffymetrixNetAffxCsvFile() to AromaUgpFile.
# 2007-03-04
# o Added findByChipType() and fromChipType().
# o Now the default path for createFromCdf() is the same as for the CDF.
# 2007-03-03
# o Now inherits from generic AromaGenomePositionFile.
# 2007-03-02
# o Created. Can import genome information data.
############################################################################
