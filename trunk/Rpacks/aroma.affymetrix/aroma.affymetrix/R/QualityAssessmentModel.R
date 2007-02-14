###########################################################################/**
# @RdocClass QualityAssessmentModel
#
# @title "The QualityAssessmentModel class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{plm}{A @see "ProbeLevelModel".}
#   \item{tags}{A @character @vector of tags.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setConstructorS3("QualityAssessmentModel", function(plm=NULL, tags="*", ...) {
  # Argument 'plm':
  if (!is.null(plm)) {
    if (!inherits(plm, "ProbeLevelModel"))
      throw("Argument 'plm' is not a ProbeLevelModel: ", class(plm)[1]);
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "QC";
  }


  extend(Object(), "QualityAssessmentModel",
    .plm = plm,
    .tags = tags
  )
})


setMethodS3("as.character", "QualityAssessmentModel", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, "Chip-effect set:");
  s <- c(s, paste("   ", as.character(getChipEffects(this))));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getDataSet", "QualityAssessmentModel", function(this, ...) {
  getDataSet(getPlm(this), ...);
})

setMethodS3("getPlm", "QualityAssessmentModel", function(this, ...) {
  this$.plm;
})

setMethodS3("getChipEffects", "QualityAssessmentModel", function(this, ...) {
  getChipEffects(getPlm(this));
})

setMethodS3("nbrOfArrays", "QualityAssessmentModel", function(this, ...) {
  ces <- getChipEffects(this);
  nbrOfArrays(ces);
})

setMethodS3("getCdf", "QualityAssessmentModel", function(this, ...) {
  ces <- getChipEffects(this);
  getCdf(ces);
}, protected=TRUE)

setMethodS3("getName", "QualityAssessmentModel", function(this, ...) {
  getName(this$.plm);
})

setMethodS3("getTags", "QualityAssessmentModel", function(this, ...) {
  ces <- getChipEffects(this);
  c(getTags(ces), this$.tags);
})

setMethodS3("getFullName", "QualityAssessmentModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})

setMethodS3("getRootPath", "QualityAssessmentModel", function(this, ...) {
  "qcData";
})

setMethodS3("getPath", "QualityAssessmentModel", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  cdf <- getCdf(this);
  chipType <- getChipType(cdf, fullname=FALSE);

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
# @RdocMethod getResiduals
#
# @title "Calculates the residuals from a probe-level model"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{}{}
# }
#
# \value{
#   Returns an @see "QualityAssessmentSet".
# }
#
# \details{
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getResiduals", "QualityAssessmentModel", function(this, units=NULL, path=NULL, name="qcData", tags="*", ram=1, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  resFcn <- function(kk) {
    yL <- .subset2(rawDataList, kk);
    thetaL <- .subset2(chipEffectList, kk);
    phiL <- .subset2(probeAffinityList, kk);
    nbrOfGroups <- length(yL);
    res <- lapply(seq_len(nbrOfGroups), FUN=function(gg) {
      y <- .subset2(.subset2(yL, gg), "intensities");
      theta <- .subset2(.subset2(thetaL, gg), "theta")[1,];
      phi <- .subset2(.subset2(phiL, gg), "phi");
      yhat <- outer(phi, theta, FUN="*");
      eps <- (y - yhat);
      list(eps=eps);
    })
    names(res) <- names(yL);
    res;
  } # resFcn()
    

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (is.null(path)) {
    cdf <- getCdf(this);
    chipType <- getChipType(cdf, fullname=FALSE);
    path <- file.path(name, getName(this), "residuals", chipType);
  }
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path);
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "residuals";
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  
  unitsPerChunk <- ram * 100000/nbrOfArrays(this);
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  # If residuals already calculated, and if force==FALSE, just return
  # a CelSet with the previous calculations

  verbose && enter(verbose, "Calculating PLM residuals");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get data and parameter objects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  plm <- getPlm(this);
  ds <- getDataSet(plm);
  ces <- getChipEffects(plm);
  paf <- getProbeAffinities(plm);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathnames
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mkdirs(path);
  filenames <- getNames(ces);
  pathnames <- sapply(filenames, function(filename) {
    filename <- sprintf("%s,residuals.cel", filename);
    Arguments$getWritablePathname(filename, path=path);
  });
  nbrOfFiles <- length(pathnames);
  

  nbrOfArrays <- nbrOfArrays(ds);
  cdf <- getCdf(ds);
  cdfHeader <- getHeader(cdf);  # Might be used later

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # find units to do; first check whether last file in list exists
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (isFile(pathnames[nbrOfFiles])) {
    qaf <- QualityAssessmentFile$fromFile(pathnames[nbrOfFiles]);
    unitsToDo <- findUnitsTodo(qaf, units=units, ..., verbose=less(verbose));
  } else {
    if (is.null(units)) {
      unitsToDo <- units;
    } else {
      unitsToDo <- seq(length=nbrOfUnits);
    }
  }

  nbrOfUnits <- length(unitsToDo);
  verbose && printf(verbose, "Number of units to do: %d\n", nbrOfUnits);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate residuals in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                    nbrOfChunks, unitsPerChunk);
  head <- seq(length=unitsPerChunk);
  count <- 1;
  while (length(unitsToDo) > 0) {
    if (length(unitsToDo) < unitsPerChunk) {
      head <- 1:length(unitsToDo);
    }
    units <- unitsToDo[head];
    verbose && printf(verbose, "Chunk #%d of %d (%d units)\n",
                                        count, nbrOfChunks, length(units));

    verbose && enter(verbose, "Getting cell indices");
    cdfList <- getCellIndices(cdf, units=units, stratifyBy="pm", ...);
    verbose && str(verbose, cdfList[[1]]);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Retrieving raw data");
    rawDataList <- readUnits(ds, units=units, stratifyBy="pm", verbose=less(verbose));
    verbose && str(verbose, rawDataList[[1]]);
    verbose && exit(verbose);

    verbose && enter(verbose, "Retrieving chip-effect estimates");
    chipEffectList <- readUnits(ces, units=units, verbose=less(verbose));
    verbose && str(verbose, chipEffectList[[1]]);
    verbose && exit(verbose);

    verbose && enter(verbose, "Retrieving probe-affinity estimates");
    probeAffinityList <- readUnits(paf, units=units, verbose=less(verbose));
    verbose && str(verbose, probeAffinityList[[1]]);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating residuals");
    residualsList <- lapply(head, FUN=resFcn);
    verbose && str(verbose, residualsList[[1]]);
    verbose && exit(verbose);
    
    # Store residuals
    for (kk in seq(along=pathnames)) {
      # Back-transform data to intensity scale and encode as CEL structure
      verbose && enter(verbose, "Encode as CEL structure");
      data <- lapply(residualsList, FUN=function(groups) {
        lapply(groups, FUN=function(group) {
          eps <- .subset2(group, "eps")[,kk];
          ones <- rep(1, length=length(eps));
          list(intensities=eps, stdvs=ones, pixels=ones);
        })
      })
      verbose && str(verbose, data[[1]]);
      verbose && exit(verbose);

      pathname <- pathnames[kk];
      if (!isFile(pathname)) {
        df <- getFile(ds, kk);
        celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=getName(df));
        createCel(pathname, header=celHeader, verbose=less(verbose));
      }

      verbose && enter(verbose, "updating file #", kk);
      updateCelUnits(pathname, cdf=cdfList, data=data);
      verbose && exit(verbose);
    } # for (kk ...)
    
    unitsToDo <- unitsToDo[-head];
    count <- count + 1;
  } # while (length(unitsToDo) > 0)

  # Return residual set
  res <- QualityAssessmentSet$fromFiles(path=path, pattern=",residuals.[cC][eE][lL]$");

  verbose && exit(verbose);

  res;
})



###########################################################################/**
# @RdocMethod getWeights
#
# @title "Calculates the weights from the robust fit to a probe-level model"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{}{}
# }
#
# \value{
#   Returns an @see "QualityAssessmentSet".
# }
#
# \details{
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getWeights", "QualityAssessmentModel", function(this, path=NULL, name="qcData", tags="*", ram=1, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  resFcn <- function(kk) {
    yL <- .subset2(rawDataList, kk);
    thetaL <- .subset2(chipEffectList, kk);
    phiL <- .subset2(probeAffinityList, kk);
    nbrOfGroups <- length(yL);
    res <- lapply(nbrOfGroups, FUN=function(gg) {
      y <- .subset2(.subset2(yL, gg), "intensities");
      theta <- .subset2(.subset2(thetaL, gg), "theta")[1,];
      phi <- .subset2(.subset2(phiL, gg), "phi");
      yhat <- outer(phi, theta, FUN="+");
      eps <- (y - yhat);
      mad <- 1.4826 * median(abs(yhat));
      matrix(psi.huber(yhat/mad), ncol=ncol(y));
    })
    res;
  } # resFcn()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (is.null(path)) {
    cdf <- getCdf(this);
    chipType <- getChipType(cdf, fullname=FALSE);
    path <- file.path(name, getName(this), "weights", chipType);
  }
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path);
  }

  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "weights";
  }

  # Argument 'unitsPerChunk':
  unitsPerChunk <- ram * 100000/nbrOfArrays(this);
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  # If residuals already calculated, and if force==FALSE, just return
  # a CelSet with the previous calculations

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  names <- getNames(getChipEffects(this));
  pathname <- sapply(names, function(name) {
    filename <- sprintf("%s,weights.cel", name);
    Arguments$getWritablePathname(filename, path=path)
  });
  nbrOfFiles <- length(pathname);
  
  ds <- getDataSet(this);
  ces <- getChipEffects(this);
  paf <- getProbeAffinities(getPlm(this));
  
  nbrOfUnits <- nbrOfUnits(getCdf(ces));

  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

# find number of units to do; first check whether last file in list
# exists

  if (isFile(pathname[nbrOfFiles])) {
    qaf <- QualityAssessmentFile$fromFile(pathname[nbrOfFiles]);
    unitsToDo <- findUnitsTodo(qaf);
  } else {
    unitsToDo <- 1:nbrOfUnits;
  }

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                    nbrOfChunks, unitsPerChunk);
  head <- 1:unitsPerChunk;
  verbose && enter(verbose, "Extracting unit data");
  count <- 1;
  while (length(unitsToDo) > 0) {
    if (length(unitsToDo) < unitsPerChunk) {
      head <- 1:length(unitsToDo);
    }
    units <- unitsToDo[head];
    verbose && printf(verbose, "Chunk #%d of %d (%d units)\n",
                                        count, nbrOfChunks, length(units));

    rawDataList <- readUnits(ds, units=units, transforms=list(log2), verbose=less(verbose), stratifyBy="pm");
    chipEffectList <- readUnits(ces, units=units, transforms=list(log2), verbose=less(verbose));
    probeAffinityList <- readUnits(paf, units=units, transforms=list(log2), verbose=verbose);

    
   weightsList <- lapply(head, FUN=resFcn);

# update output files
    
    cdf <- getCellIndices(getCdf(ds), units=units, stratifyBy="pm", ...);
  
    for (kk in seq(pathname)) {
      if (!isFile(pathname[kk])) {
        cdfHeader <- getHeader(getCdf(ds));
        celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=getName(getFile(ds,kk)));
        createCel(pathname[kk], header=celHeader, verbose=less(verbose));
      }
      data <- lapply(weightsList, function(x){
        nrow <- nrow(x); 
        list(list(
          intensities=2^x[,kk], 
          stdvs=rep(1, nrow), 
          pixels=rep(1, nrow)
        ))
      });
      updateCelUnits(pathname[kk], cdf=cdf, data=data);
    }
    
    unitsToDo <- unitsToDo[-head];
    count <- count + 1;
  }

  res <- QualityAssessmentSet$fromFiles(path=path, pattern=",weights.[cC][eE][lL]$");
  setAlias(res, getName(this));

  res;
})




setMethodS3("plotNuse", "QualityAssessmentModel", function(this, subset=NULL, verbose=FALSE, main="NUSE", ...) {

# ... : additional arguments to bxp().

  verbose <- Arguments$getVerbose(verbose);
  
  plm <- getPlm(this);
  ces <- getChipEffects(plm);
  cdfMono <- getCdf(ces);
  nbrOfUnits <- nbrOfUnits(cdfMono);
  
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, range=c(1, nbrOfUnits));
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }

  # get the vector of median stdvs
  verbose && enter(verbose, "Extracting standard errors");
  avg <- getAverageLog(getChipEffects(plm), field="stdvs", indices=units, verbose=verbose)
  verbose && exit(verbose);
  medianSE <- getData(avg, indices=units, "intensities")$intensities;
  medianSE <- log2(medianSE);

  verbose && enter(verbose, "Calculating summaries for ", nbrOfArrays(this), " arrays");
  
  # for each file, sweep through and calculate statistics for boxplot.
  boxplotStats <- list();
  for (kk in seq(ces)) {
    stdvs <- getData(getFile(ces, kk), indices=units, "stdvs")$stdvs;
    stdvs <- log2(stdvs);
    boxplotStats[[kk]] <- boxplot.stats(stdvs/medianSE);
  }
  rm(avg);
  rm(stdvs);    # not needed any more

  verbose && exit(verbose);
  
  # make a new list from boxplotStats which has correct structure to
  # pass to bxp().

  bxpStats <- list();

  elementNames <- names(boxplotStats[[1]]);

  for (name in elementNames) {
    suppressWarnings(bxpStats[[name]] <- do.call("cbind", lapply(boxplotStats, function(x){x[[name]]})));
  }
  
  bxp(bxpStats, main=main, ...);
})


# ... : additional arguments to bxp().
setMethodS3("plotRle", "QualityAssessmentModel", function(this, subset=NULL, verbose=FALSE, main="RLE", ...) {
  verbose <- Arguments$getVerbose(verbose);

  plm <- getPlm(this);
  ces <- getChipEffects(plm);
  cdfMono <- getCdf(ces);
  nbrOfUnits <- nbrOfUnits(cdfMono);
  
  # Argument 'subset':
  if (!(is.null(subset))) {
    getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1);
    if (!getFraction) {
      units <- Arguments$getIndices(subset, range=c(1, nbrOfUnits));
    } else {
      units <- identifyCells(cdfMono, indices=subset, verbose=less(verbose));
    }
  } else {
    units <- 1:nbrOfUnits;
  }

  # get the vector of median stdvs
  verbose && enter(verbose, "Extracting chip effects");
  avg <- getAverageLog(ces, field="intensities", indices=units, verbose=verbose)
  verbose && exit(verbose);
  medianLE <- getData(avg, indices=units, "intensities")$intensities;
  medianLE <- log2(medianLE);

  verbose && enter(verbose, "Calculating summaries for ", nbrOfArrays(this), " arrays");
  
  # for each file, sweep through and calculate statistics for boxplot.
  boxplotStats <- list();
  for (kk in seq(ces)) {
    theta <- getData(getFile(ces, kk), indices=units, "intensities")$intensities;
    theta <- log2(theta);
    boxplotStats[[kk]] <- boxplot.stats(theta-medianLE);
  }
  rm(avg, theta, medianLE, units);  # not needed any more

  verbose && exit(verbose);

  
  # make a new list from boxplotStats which has correct structure to
  # pass to bxp().
  bxpStats <- list();
  elementNames <- names(boxplotStats[[1]]);
  for (name in elementNames) {
    suppressWarnings({
      args <- lapply(boxplotStats, FUN=function(x) { x[[name]] });
      bxpStats[[name]] <- do.call("cbind", args);
    });
  }
  
  bxp(bxpStats, main=main, ...);
})


##########################################################################
# HISTORY:
# 2007-01-31 /HB
# o Removed any argument 'unitsPerChunk' and renamed 'moreUnits' to 'ram'.
# o Replaced the usage of getDataSet() with other alternatives, in order
#   to minimize the requirement for knowing the data set, e.g. we could
#   reload much of the PLM estimates from files.
# 2007-01-15
# o Renamed to QualityAssessmentModel from QcInfo.
# 2007-01-13
# o Added getWeights().
# 2007-01-06
# o Created.
##########################################################################
