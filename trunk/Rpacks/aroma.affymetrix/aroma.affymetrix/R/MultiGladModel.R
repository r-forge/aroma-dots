###########################################################################/**
# @RdocClass MultiGladModel
#
# @title "The MultiGladModel class"
#
# \description{
#  @classhierarchy
#
#  This class is similar to the @see "GladModel" class, but it can handle
#  multiple chip types.  For now, consider this class to be very much
#  a trial-and-error class.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesList}{A @list of @see "CnChipEffectSet" objects.}
#   \item{referenceList}{A @see "CnChipEffectFile" objects.}
#   \item{tags}{A @character @vector of tags.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#  @see "GladModel".
# }
#*/###########################################################################
setConstructorS3("MultiGladModel", function(cesList=NULL, referenceList=NULL, tags="*", ...) {
  require(GLAD) || throw("Package not loaded: GLAD");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cesList':
  if (!is.null(cesList)) {
    if (!is.list(cesList)) {
      cesList <- list(cesList);
    }

    for (ces in cesList) {
      if (!inherits(ces, "CnChipEffectSet")) {
        throw("Argument 'cesList' contains a non-CnChipEffectSet object: ", class(ces));
      }

      # Currently only total copy-number estimates are accepted
      if (!ces$combineAlleles) {
        throw("Unsupported copy-number chip effects. Currently only total copy-number estimates are supported: ces$combineAlleles == FALSE");
      }
    }
  }

  # Argument 'referenceList':
  if (!is.null(referenceList)) {
    if (!is.list(referenceList)) {
      referenceList <- list(referenceList);
    }

    if (length(referenceList) != length(cesList)) {
      throw("The number of reference files does not match the number of chip-effect sets: ", length(cesList), " != ", length(referenceList));
    }

    # Validate consistency between the chip-effect sets and the reference files
    for (kk in seq(along=cesList)) {
      ces <- cesList[[kk]];
      ref <- referenceList[[kk]];

      if (!inherits(ref, "CnChipEffects")) {
        throw("Argument 'referenceList' contains a non-CnChipEffects object: ",
                                                            class(ref));
      }

      if (ref$combineAlleles != ces$combineAlleles) {
         throw("The reference chip effects are not compatible with the chip-effect set. One is combining the alleles the other is not.");
      }

      if (ref$mergeStrands != ces$mergeStrands) {
         throw("The reference chip effects are not compatible with the chip-effect set. One is merging the strands the other is not.");
      }
    } # for (kk in ...)
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "GLAD";
  }


  

  extend(Object(), "MultiGladModel",
    .cesList = cesList,
    .referenceList = referenceList
  )
})


setMethodS3("as.character", "MultiGladModel", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  nbrOfChipTypes <- nbrOfChipTypes(this);
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes));
  s <- c(s, "Chip-effect set & reference file pairs:");
  cesList <- getListOfChipEffects(this);
  refList <- getListOfReferences(this);
  for (kk in seq(length=nbrOfChipTypes)) {
    s <- c(s, sprintf("Pair #%d:", kk));
    s <- c(s, "Chip-effect set:");
    ces <- cesList[[kk]];
    ref <- refList[[kk]];
    s <- c(s, as.character(ces));
    s <- c(s, "Reference file:");
    if (is.null(ref)) {
      s <- c(s, "<average across arrays>");
    } else {
      s <- c(s, as.character(ref));
    }
  }
#  s <- c(s, "Genome information:", as.character(getGenomeInformation(this)));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})

setMethodS3("getListOfChipEffects", "MultiGladModel", function(this, ...) {
  this$.cesList;
})


setMethodS3("nbrOfChipTypes", "MultiGladModel", function(this, ...) {
  length(getListOfChipEffects(this, ...));
})

setMethodS3("getListOfReferences", "MultiGladModel", function(this, ...) {
  res <- this$.referenceList;
  if (is.null(res)) {
    res <- vector("list", nbrOfChipTypes(this));
  }
  res;
})

setMethodS3("setListOfReferences", "MultiGladModel", function(this, referenceList=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'referenceList':
  if (!is.null(referenceList)) {
    if (!is.list(referenceList)) {
      referenceList <- list(referenceList);
    }

    cesList <- getListOfChipEffects(this);
    if (length(referenceList) != length(cesList)) {
      throw("The number of reference files does not match the number of chip-effect sets: ", length(cesList), " != ", length(referenceList));
    }

    # Validate consistency between the chip-effect sets and the reference files
    for (kk in seq(along=cesList)) {
      ces <- cesList[[kk]];
      ref <- referenceList[[kk]];

      if (!inherits(ref, "CnChipEffects")) {
        throw("Argument 'referenceList' contains a non-CnChipEffects object: ",
                                                            class(ref));
      }

      if (ref$combineAlleles != ces$combineAlleles) {
         throw("The reference chip effects are not compatible with the chip-effect set. One is combining the alleles the other is not.");
      }

      if (ref$mergeStrands != ces$mergeStrands) {
         throw("The reference chip effects are not compatible with the chip-effect set. One is merging the strands the other is not.");
      }
    } # for (kk in ...)
  }

  this$.referenceList <- referenceList;
})


setMethodS3("getListOfCdfs", "MultiGladModel", function(this, ...) {
  cesList <- getListOfChipEffects(this);
  lapply(cesList, FUN=getCdf);
}, protected=TRUE)


setMethodS3("getChipTypes", "MultiGladModel", function(this, merge=FALSE, collapse="+", ...) {
  cdfList <- getListOfCdfs(this);
  chipTypes <- sapply(cdfList, FUN=getChipType);

  chipTypes <- sapply(chipTypes, FUN=function(s) {
    gsub("[,-]monocell", "", s);
  })

  if (merge) {
    chipTypes <- mergeByCommonTails(chipTypes, collapse=collapse);
  }

  chipTypes;
}, protected=TRUE)


setMethodS3("getChipType", "MultiGladModel", function(this, ...) {
  getChipTypes(this, merge=TRUE, ...);
})


setMethodS3("getTableOfArrays", "MultiGladModel", function(this, ...) {
  cesList <- getListOfChipEffects(this);

  # Get all chip types for this data set
  chipTypes <- getChipTypes(this);
  nbrOfChipTypes <- length(cesList);

  # Get array names
  arrayNames <- lapply(cesList, FUN=getNames);
  names(arrayNames) <- chipTypes;

  # Get all unique array names
  allNames <- unlist(arrayNames, use.names=FALSE);
  allNames <- unique(allNames);

  # Create table of arrays
  nbrOfArrays <- length(allNames);
  X <- matrix(NA, nrow=nbrOfArrays, ncol=nbrOfChipTypes);
  dimnames(X) <- list(allNames, chipTypes);
  for (chipType in chipTypes) {
    names <- arrayNames[[chipType]];
    idx <- match(names, allNames);
    X[idx,chipType] <- seq(along=idx);
  }

  X;
})

setMethodS3("getNames", "MultiGladModel", function(this, ...) {
  rownames(getTableOfArrays(this, ...));
})

setMethodS3("getArrays", "MultiGladModel", function(this, ...) {
  getNames(this, ...);
})

setMethodS3("nbrOfArrays", "MultiGladModel", function(this, ...) {
  length(getNames(this, ...));
})


setMethodS3("getName", "MultiGladModel", function(this, collapse="+", ...) {
  # Get name of chip-effect sets
  cesList <- getListOfChipEffects(this);

  # Get names
  names <- lapply(cesList, FUN=getName);
  names <- unlist(names, use.names=FALSE);

  # Merge names
  names <- mergeByCommonTails(names, collapse=collapse);

  names;
})


setMethodS3("getReferenceName", "MultiGladModel", function(this, collapse="+", ...) {
  # Get name of the reference files
  refList <- getListOfReferences(this);

  # Get names
  names <- lapply(refList, FUN=function(file) {
    if (is.null(file)) NULL else getName(file);
  });
  names <- unlist(names, use.names=FALSE);

  # Merge names
  names <- mergeByCommonTails(names, collapse=collapse);

  names;
})


setMethodS3("getTags", "MultiGladModel", function(this, ...) {
  # Get tags of chip-effect set
  cesList <- getListOfChipEffects(this);

  # Get data set tags
  tags <- lapply(cesList, FUN=getTags);

  # Keep unique tags
  tags <- unlist(tags, use.names=FALSE);
  tags <- unique(tags);

  # Add model tags
  c(tags, this$.tags);
})


setMethodS3("getFullName", "MultiGladModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "MultiGladModel", function(this, ...) {
  "gladData";
})

setMethodS3("getPath", "MultiGladModel", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  chipType <- getChipType(this);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})

###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the chromosomes available"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChromosomes", "MultiGladModel", function(static, ...) {
  c(1:22, "X");
}, static=TRUE)



setMethodS3("getListOfGenomeInformations", "MultiGladModel", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving genome informations");
  cdfList <- getListOfCdfs(this);
  giList <- lapply(cdfList, getGenomeInformation, verbose=less(verbose));
  verbose && exit(verbose);

  giList;
})


setMethodS3("getChipEffectFiles", "MultiGladModel", function(this, array, ...) {
  cesList <- getListOfChipEffects(this);
  arrayTable <- getTableOfArrays(this);
  chipTypes <- colnames(arrayTable);
  idxs <- arrayTable[array,];
  res <- vector("list", length(chipTypes));
  for (kk in seq(along=cesList)) {
    idx <- idxs[kk];
    if (is.na(idx))
      next;
    ces <- cesList[[kk]];
    res[[kk]] <- getFile(ces, idx);
  }
  names(res) <- chipTypes;
  res;
})

setMethodS3("getReferenceFiles", "MultiGladModel", function(this, ...) {
  cesList <- getListOfChipEffects(this);
  refList <- getListOfReferences(this);
  for (kk in seq(length=nbrOfChipTypes(this))) {
    ces <- cesList[[kk]];
    ref <- refList[[kk]];
    if (is.null(ref)) {
      verbose && enter(verbose, "No reference specified. Calculating average chip effects");
      ref <- getAverageFile(ces, verbose=less(verbose));
      verbose && exit(verbose);
      refList[[kk]] <- ref;
    }
  }
  refList;
})



setMethodS3("fit", "MultiGladModel", function(this, arrays=1:nbrOfArrays(this), chromosomes=getChromosomes(this), ..., .retResults=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- Arguments$getIndices(arrays, range=c(1,nbrOfArrays(this)));

  # Argument 'chromosomes':
  chromosomes <- Arguments$getCharacters(chromosomes);
#  chromosomes[chromosomes == "23"] <- "X";
  chromosomes <- intersect(chromosomes, getChromosomes(this));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);
  mkdirs(path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get chip effects
  cesList <- getListOfChipEffects(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving reference chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  refList <- getReferenceFiles(this);
  verbose && cat(verbose, "Using references:");
  verbose && print(verbose, refList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get reference annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get reference tags across chip types
  tags <- lapply(refList, FUN=function(file) {
    if (is.null(file)) NULL else getTags(file);
  });
  tags <- unlist(tags, use.names=FALSE);
  tags <- unique(tags);
  # Exclude unwanted tags
  refTags <- setdiff(tags, "chipEffects");

  # Add combined reference name
  refName <- getReferenceName(this);
  refTags <- c(refName, refTags);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  arrayNames <- getArrays(this)[arrays];
  nbrOfArrays <- length(arrayNames);
  for (aa in seq(length=nbrOfArrays)) {
    array <- arrays[aa];
    arrayName <- arrayNames[aa];

    ceList <- getChipEffectFiles(this, array=array);

    # Get chip-effect tags across chip types
    tags <- lapply(ceList, FUN=function(file) {
      if (is.null(file)) NULL else getTags(file);
    });
    tags <- unlist(tags, use.names=FALSE);
    tags <- unique(tags);
    # Exclude unwanted tags
    ceTags <- setdiff(tags, "chipEffects");


    res[[arrayName]] <- list();
    for (chr in chromosomes) {
      chrIdx <- match(chr, c(1:22, "X", "Y"));
      verbose && enter(verbose, 
                             sprintf("Array %s (#%d of %d) on chromosome %s", 
                                           arrayName, aa, nbrOfArrays, chr));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get pathname
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Add tags chrNN,glad,<reference tags>
      tags <- c(ceTags, sprintf("chr%02d", chrIdx));
      tags <- c(tags, "glad", refTags);
      fullname <- paste(c(arrayName, tags), collapse=",");
      filename <- sprintf("%s.xdr", fullname);
      pathname <- filePath(path, filename);

      # Already done?
      if (isFile(pathname)) {
        verbose && enter(verbose, "Loading results from file");
        verbose && cat(verbose, "Pathname: ", pathname);
        fit <- loadObject(pathname);
        verbose && exit(verbose);
      } else {
        fit <- fitOne(this, ceList=ceList, refList=refList, chromosome=chr, 
                                                 ..., verbose=less(verbose));
        verbose && enter(verbose, "Saving to file");
        verbose && cat(verbose, "Pathname: ", pathname);
        saveObject(fit, file=pathname);
        verbose && exit(verbose);
      }

      verbose && enter(verbose, "Calling onFit() hooks");
      callHooks("onFit.fitGlad.MultiGladModel", fit=fit, fullname=fullname);
      verbose && exit(verbose);

      if (.retResults)
        res[[arrayName]][[chr]] <- fit;

      rm(fit);

      verbose && exit(verbose);
    } # for (chr in ...)
  } # for (aa in ...)

  invisible(res);
})



setMethodS3("plot", "MultiGladModel", function(x, ..., pixelsPerMb=3, zooms=2^(0:7), pixelsPerTick=2.5, height=400, imageFormat="current", skip=TRUE, path=NULL, callList=NULL, verbose=FALSE) {
  # To please R CMD check.
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pixelsPerMb':
  pixelsPerMb <- Arguments$getDouble(pixelsPerMb, range=c(0.001,9999));

  # Argument 'zooms':
  zooms <- Arguments$getIntegers(zooms, range=c(1,9999));
  zooms <- unique(zooms);

  # Argument 'pixelsPerMb':
  pixelsPerTick <- Arguments$getDouble(pixelsPerTick, range=c(1,256));

  # Argument 'height':
  height <- Arguments$getInteger(height, range=c(1,4096));

  # Argument 'callList':
  chipTypes <- getChipTypes(this);
  if (length(callList) > 0) {
    if (!is.list(callList))
      callList <- list(callList);

    if (length(callList) != length(chipTypes)) {
      throw("Number of elements in argument 'callList' does not match the number of chip types: ", length(callList), " != ", length(chipTypes));
    }

    if (is.null(names(callList)))
      names(callList) <- chipTypes;

    for (chipType in chipTypes) {
      callSet <- callList[[chipType]];
      if (!is.null(callSet)) {
        if (!inherits(callSet, "GenotypeCallSet"))
          throw("Argument 'callList' contains a non-GenotypeCallSet: ", class(callSet)[1]);

        if (getChipType(callSet) != chipType) {
          throw("Argument 'callList' contains a GenotypeCallSet for a different chip type than the corresponding chip-effect set: ", getChipType(callSet), " != ", chipType);
        }
      }
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get genome annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  genome <- readTable("annotations/hgChromosomes.txt", header=TRUE, 
                            colClasses=c(nbrOfBases="integer"), row.names=1);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Output path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rootPath <- "chromosomeExplorer";

  # Get chip type
  chipType <- getChipType(this);

  # The figure path
  if (is.null(path)) {
    modelName <- paste(getTags(this), collapse=",");
    path <- filePath(rootPath, getName(this), modelName, chipType, expandLinks="any");
    path <- filePath(path, expandLinks="any");
  }
  mkdirs(path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the PNG device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(imageFormat)) {
    imageFormat <- "current";
  }

  if (identical(imageFormat, "current")) {
    plotDev <- NULL;
    zooms <- zooms[1];
  } else if (identical(imageFormat, "screen")) {
    screenDev <- function(pathname, width, height, ..., xpinch=50, ypinch=xpinch) {
      # Dimensions are in pixels. Rescale to inches
      width <- width/xpinch;
      height <- height/ypinch;
      x11(width=width, height=height, ...);
    }

    # When plotting to the screen, use only the first zoom
    zooms <- zooms[1];
    plotDev <- screenDev;
  } else if (identical(imageFormat, "png")) {
    pngDev <- System$findGraphicsDevice();
    plotDev <- pngDev;
  }


  callCols <- c("-"="lightgray", AA="red", AB="blue", BB="red", NC="orange");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Define the plot function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  on.exit({
    setHook("onFit.fitGlad.MultiGladModel", NULL, action="replace");
  })

  setHook("onFit.fitGlad.MultiGladModel", function(fit, fullname) {
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    tryCatch({
    # Extract the array name from the full name
    arrayName <- gsub("^([^,]*).*$", "\\1", fullname);

    # Extract the chromsome from the GLAD fit object
    pv <- fit$profileValues;
    chromosome <- unique(pv$Chromosome);  # Should only be one!
    chromosome[chromosome == "23"] <- "X";

    # Infer the length (in bases) of the chromosome
    chromosomeIdx <- match(chromosome, getChromosomes(this));
    nbrOfBases <- genome$nbrOfBases[chromosomeIdx];
    widthMb <- nbrOfBases / 1e6;

    verbose && enter(verbose, sprintf("Plotting %s for chromosome %s (%02d) [%.2fMb]", arrayName, chromosome, chromosomeIdx, widthMb));

    for (zz in seq(along=zooms)) {
      zoom <- zooms[zz];

      # Create the pathname to the file
      imgName <- sprintf("%s,glad,chr%02d,x%04d.%s", arrayName, chromosomeIdx, zoom, imageFormat);
      pathname <- filePath(path, imgName);

      # pngDev() (that is bitmap()) does not accept spaces in pathnames
      pathname <- gsub(" ", "_", pathname);
      if (!imageFormat %in% c("screen", "current")) {
        if (skip && isFile(pathname)) {
          next;
        }
      }
      # Calculate Mbs per ticks
      ticksBy <- 10^ceiling(log10(pixelsPerTick / (zoom * pixelsPerMb)));
      width <- as.integer(zoom * widthMb * pixelsPerMb);
      # Plot to PNG file
      verbose && printf(verbose, "Pathname: %s\n", pathname);
      verbose && printf(verbose, "Dimensions: %dx%d\n", width, height);
      verbose && printf(verbose, "Ticks by: %f\n", ticksBy);

      if (!is.null(plotDev))
        plotDev(pathname, width=width, height=height);

      tryCatch({
        verbose && enter(verbose, "Plotting graph");
        opar <- par(xaxs="r");
        plot(fit, ticksBy=ticksBy);
        stext(chipType, side=4, pos=1, line=0, cex=0.7, col="gray");

        if (!is.null(callList)) {
          verbose && enter(verbose, "Adding genotype calls");

          # Get (chip type, unit) information
          chipType <- pv$chipType;
          unit <- pv$unit;

          # Figure out where to put the genotype track
          ylim <- par("usr")[3:4];
          ylim <- ylim + c(+1,-1)*0.04*diff(ylim);
          ylim <- ylim + c(+1,-1)*0.04*diff(ylim);

          for (chipType in chipTypes) {
            # Got genotype calls for this chip type?
            callSet <- callList[[chipType]];
            if (is.null(callSet))
              next;

            # Got copy-number estimates for this chip type?
            idxs <- which(pv$chipType == chipType);
            if (length(idxs) == 0)
              next;

            # Got copy-number estimates for this array?
            if (!arrayName %in% getNames(callSet))
              next;

            # Get the subset of genotype calls for this array and chromosome.
            units <- pv$unit[idxs];
            call <- callSet[units, arrayName];
            call <- as.character(call);
            # Get the positions of these calls
            x <- pv$PosBase[idxs]; 
            x <- x/1e6;

            # Plot the homozygote/heterozygote genotype tracks
            y <- rep(ylim[1], length(callCols));
            names(y) <- names(callCols);
            y["AB"] <- y["AB"] + 0.02*diff(ylim);
            y <- y[call];
            points(x,y, pch=".", cex=2, col=callCols[call]);
 
            rm(idxs, call, callSet, units, x, y);  # Not needed anymore
          } # for (chipType ...)
          verbose && exit(verbose);
        } # if (!is.null(callList))

        verbose && exit(verbose);
      }, error = function(ex) {
        print(ex);
      }, finally = {
        par(opar);
        if (!imageFormat %in% c("screen", "current"))
          dev.off();
      });
    } # for (zoom ...)

    verbose && exit(verbose);
    }, error = function(ex) {
      cat("ERROR caught in onFit.fitGlad.MultiGladModel():\n");
      print(ex);
    }) # tryCatch()
  }, action="replace")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Start fitting and plotting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit(this, ..., .retResults=FALSE, verbose=verbose);

  invisible();
})


setMethodS3("getLog2Ratios", "MultiGladModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the regions for each of the GLAD fits (per array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  suppressWarnings({
    res <- fit(this, ..., verbose=less(verbose,10));
  })

  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    for (fit in arrayFits) {
      suppressWarnings({
        df0 <- getRegions(fit, ...);
      })
      df <- rbind(df, df0);
    }
    rownames(df) <- seq(length=nrow(df));
  })

  res;
}, protected=TRUE) # getLog2Ratios()



setMethodS3("getRegions", "MultiGladModel", function(this, ..., url="ucsc", margin=10, flat=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'url':
  if (identical(url, "ucsc")) {
    url <- function(chromosome, start, stop) {
      uri <- "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=vertebrate&org=Human&db=hg18";
      sprintf("%s&position=chr%s%%3A%d-%d", uri, chromosome, as.integer(start), as.integer(stop));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the regions for each of the GLAD fits (per array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  suppressWarnings({
    res <- fit(this, ..., verbose=less(verbose,10));
  })

  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    for (fit in arrayFits) {
      suppressWarnings({
        df0 <- getRegions(fit, ...);
      })
      df <- rbind(df, df0);
    }
    rownames(df) <- seq(length=nrow(df));

    # Add URL?
    if (!is.null(url)) {
      chromosome <- df[,"Chromosome"];
      start <- df[,"start"];
      stop <- df[,"stop"];
      m <- margin*abs(stop-start);
      start <- start-m;
      start[start < 0] <- 0;
      stop <- stop + m;
      urls <- character(nrow(df));
      for (rr in seq(along=urls)) { 
        urls[rr] <- url(chromosome[rr], start[rr], stop[rr]);
      }
      df <- cbind(df, url=urls);
    }

    df;
  })

  if (flat) {
    df <- NULL;
    for (kk in seq(along=res)) {
      df <- rbind(df, cbind(sample=names(res)[kk], res[[kk]]));
      res[[kk]] <- NA;
    }
    row.names(df) <- seq(length=nrow(df));
    res <- df;
  }  

  res;
})


setMethodS3("writeRegions", "MultiGladModel", function(this, arrays=1:nbrOfArrays(this), format=c("xls", "wig"), digits=3, ..., oneFile=TRUE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- Arguments$getIndices(arrays, range=c(1,nbrOfArrays(this)));

  # Argument 'format':
  format <- match.arg(format);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Setup
  fullname <- getFullName(this);
  arrayNames <- getArrays(this);

  path <- getPath(this);
  mkdirs(path);

  if (oneFile) {
    filename <- sprintf("%s,GLAD,regions.%s", fullname, format); 
    pathname <- filePath(path, filename);
    pathname <- Arguments$getWritablePathname(pathname);
    if (!skip && isFile(pathname))
      file.remove(pathname);
  }

  res <- list();
  for (array in arrays) {
    name <- arrayNames[array];
    verbose && enter(verbose, sprintf("Array #%d (of %d) - %s", 
                                             array, length(arrays), name));
    df <- getRegions(this, arrays=array, ..., verbose=less(verbose))[[1]];
    names(df) <- gsub("Smoothing", "log2CN", names(df));

    if (identical(format, "xls")) {
      # Append column with sample names
      df <- cbind(sample=name, df);
    } else if (identical(format, "wig")) {
      # Write a four column WIG/BED table
      df <- df[,c("Chromosome", "start", "stop", "log2CN")];

      # In the UCSC Genome Browser, the maximum length of one element
      # is 10,000,000 bases.  Chop up long regions in shorter contigs.
      verbose && enter(verbose, sprintf("Chopping up too long segment"));
      MAX.LENGTH = 10e6-1;
      start <- df[,"start"];
      stop <- df[,"stop"];
      len <- stop-start;
      tooLong <- which(len > MAX.LENGTH);
      if (length(tooLong) > 0) {
        dfXtra <- NULL;
        for (rr in tooLong) {
          x0 <- start[rr];
          while (x0 < stop[rr]) {
            x1 <- min(x0 + MAX.LENGTH, stop[rr]);
            df1 <- df[rr,];
            df1[,"start"] <- x0;
            df1[,"stop"] <- x1;
            dfXtra <- rbind(dfXtra, df1);          
            x0 <- x1+1;
          }
        }
        df <- df[-tooLong,];
        df <- rbind(df, dfXtra);
        rm(dfXtra);
        row.names(df) <- seq(length=nrow(df));
      }
      verbose && exit(verbose);
      # Make sure the items are ordered correctly
      chrIdx <- as.integer(df[,"Chromosome"]);
      o <- order(chrIdx, df[,"start"]);
      df <- df[o,];

      # All chromosomes should have prefix 'chr'.
      chrIdx <- as.integer(df[,"Chromosome"]);
      df[chrIdx == 23,"Chromosome"] <- "X";
      df[,"Chromosome"] <- paste("chr", df[,"Chromosome"], sep="");
    }

    # Apply digits
    for (cc in seq(length=ncol(df))) {
      value <- df[,cc];
      if (is.double(value)) {
        df[,cc] <- round(value, digits=digits);
      }
    }

    if (!oneFile) {
      savename <- name;
      filename <- sprintf("%s,GLAD,regions.%s", savename, format); 
      pathname <- filePath(path, filename);
      if (!oneFile && !skip && isFile(pathname))
        file.remove(pathname);
    }

    # Writing to file
    verbose && cat(verbose, "Pathname: ", pathname);
    if (identical(format, "xls")) {
      col.names <- (array == arrays[1]);
      write.table(df, file=pathname, sep="\t", col.names=col.names, row.names=FALSE, quote=FALSE, append=oneFile);
    } else if (identical(format, "wig")) {
      # Write track control
      trackAttr <- c(type="wiggle_0");
      trackAttr <- c(trackAttr, name=sprintf("\"%s\"", name));
      trackAttr <- c(trackAttr, group="\"GLAD regions\"");
      trackAttr <- c(trackAttr, priority=array);
      trackAttr <- c(trackAttr, graphType="bar");
      trackAttr <- c(trackAttr, visibility="full");
      trackAttr <- c(trackAttr, maxHeightPixels="128:96:64");
      trackAttr <- c(trackAttr, yLineOnOff="on");
# HARD WIRED FOR NOW.  TO DO /hb 2006-11-27
col <- c("117,112,179", "231,41,138");
ylim <- c(-1,1);
      if (!is.null(col)) {
        trackAttr <- c(trackAttr, color=col[1], altColor=col[2]);
      }
      if (!is.null(ylim)) {
        trackAttr <- c(trackAttr, autoScale="off", 
              viewLimits=sprintf("%.2f:%.2f ", ylim[1], ylim[2]));
      }
      trackAttr <- paste(names(trackAttr), trackAttr, sep="=");
      trackAttr <- paste(trackAttr, collapse=" ");
      trackAttr <- paste("track ", trackAttr, "\n", sep="");
      verbose && str(verbose, trackAttr);
      cat(file=pathname, trackAttr, append=oneFile);

      # Write data
      verbose && str(verbose, df);
      write.table(df, file=pathname, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=oneFile);
    }
    verbose && exit(verbose);
    res[[array]] <- df;
  }

  invisible(res);
})



##############################################################################
# HISTORY:
# 2006-12-17
# o BUG FIX: The new fitDone() in plot() choked on chr 23 (should be 'X').
# 2006-12-15
# o This class should be considered temporary, because we might design a
#   ChipEffectSet class that can contain multiple chip types, but treated as
#   if it contained one chip type, so it can be passed to the current 
#   GladModel class.  However, such a class design will require multiple 
#   inheritance etc, which will take time to develope.
# o Created from GladModel.R.
##############################################################################
