###########################################################################/**
# @RdocClass QualityAssessmentModel
#
# @title "The QualityAssessmentModel class"
#
# \description{
#  @classhierarchy
#
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
  
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "QC";
  }

  extend(Object(), "QualityAssessmentModel",
         .plm = plm,
         .tags = tags );
})

setMethodS3("getPlm", "QualityAssessmentModel", function(this, ...) {
  this$.plm;
})

setMethodS3("getChipEffects", "QualityAssessmentModel", function(this, ...) {
  getChipEffects(this$.plm);
})

setMethodS3("getName", "QualityAssessmentModel", function(this, ...) {
  getName(this$.plm);
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

setMethodS3("nbrOfArrays", "QualityAssessmentModel", function(this, ...) {
  ces <- getChipEffects(this);
  nbrOfArrays(ces);
})

setMethodS3("getCdf", "QualityAssessmentModel", function(this, ...) {
  ces <- getChipEffects(this);
  getCdf(ces);
}, protected=TRUE)

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

setMethodS3("getDataSet", "QualityAssessmentModel", function(this, ...) {
  getDataSet(getPlm(this), ...);
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
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell$", "", chipType);

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


setMethodS3("getResiduals", "QualityAssessmentModel", function(this, path=NULL, name="qcData", tags="*", unitsPerChunk=moreUnits*100000/length(getDataSet(getPlm(this))), moreUnits=1, force=FALSE, verbose=FALSE, ...) {

  # Argument 'path':
  if (is.null(path)) {
    path <- file.path(name, getName(this), "residuals", getChipType(getCdf(getDataSet(getPlm(this)))));
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
    tags[tags == "*"] <- "residuals";
  }

  # Argument 'unitsPerChunk':
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  # If residuals already calculated, and if force==FALSE, just return
  # a CelSet with the previous calculations

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filenames <- getNames(getDataSet(this));
  pathname <- sapply(filenames, function(filename){Arguments$getWritablePathname(sprintf("%s,residuals.cel", filename), path=path)});
  nbrOfFiles <- length(pathname);
  
  ds <- getDataSet(this);
  ces <- getChipEffects(this);
  paf <- getProbeAffinities(getPlm(this));
  
  nbrOfUnits <- nbrOfUnits(getCdf(ces));
  nbrOfArrays <- nbrOfArrays(this);
  
  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

# find number of units to do; first check whether last file in list
# exists

  if (isFile(pathname[nbrOfFiles])) {
    unitsToDo <- findUnitsTodo(QualityAssessmentFile$fromFile(pathname[nbrOfFiles]), ...);
  } else {
    unitsToDo <- 1:nbrOfUnits;
  }

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                    nbrOfChunks, unitsPerChunk);
  head <- 1:unitsPerChunk;
  verbose && enter(verbose, "Extracting unit data");
  count <- 1;
  idxOffset <- as.integer(0);

#  residualsList <- list();
  
  while (length(unitsToDo) > 0) {
    if (length(unitsToDo) < unitsPerChunk) {
      head <- 1:length(unitsToDo);
    }
    units <- unitsToDo[head];
    verbose && printf(verbose, "Chunk #%d of %d (%d units)\n",
                                        count, nbrOfChunks, length(units));

    transforms <- rep(list("log2"), nbrOfArrays);
    transforms <- lapply(transforms, get);
    
    rawDataList <- readUnits(ds, units=units, transforms=transforms, verbose=less(verbose), stratifyBy="pm");
    chipEffectList <- readUnits(ces, units=units, transforms=transforms, verbose=less(verbose));
    probeAffinityList <- readUnits(paf, units=units, transforms=transforms[1], verbose=verbose);

    resFcn <- function(kk) {
      rawData <- rawDataList[[kk]][[1]]$intensities;
      ce <- chipEffectList[[kk]][[1]]$theta[1,];
      pa <- probeAffinityList[[kk]][[1]]$phi;
      fittedData <- outer(pa, ce, FUN="+");
      return(rawData - fittedData);
    }
    
    verbose && enter(verbose, "Calculating residuals");

    residualsList <- lapply(head, FUN=resFcn);
    
    verbose && exit(verbose);
    
# update output files
    
    verbose && enter(verbose, "Getting cell indices");

    cdf <- getCellIndices(getCdf(ds), units=units, stratifyBy="pm", ...);

    verbose && exit(verbose);
    
    for (kk in seq(pathname)) {
      if (!isFile(pathname[kk])) {
        cdfHeader <- getHeader(getCdf(ds));
        celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=getName(getFile(ds,kk)));
        createCel(pathname[kk], header=celHeader, verbose=less(verbose));
      }
      data <- lapply(residualsList, function(x){nrow <- nrow(x); list(list(intensities=2^x[,kk], stdvs=rep(1, nrow), pixels=rep(1, nrow)))})
      verbose && enter(verbose, "updating file #", kk);
      updateCelUnits(pathname[kk], cdf=cdf, data=data);
      verbose && exit(verbose);
    }
    
    unitsToDo <- unitsToDo[-head];
    count <- count + 1;
  }

  res <- QualityAssessmentSet$fromFiles(path=path, pattern=",residuals.[cC][eE][lL]$");
  return(res);
  
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


setMethodS3("getWeights", "QualityAssessmentModel", function(this, path=NULL, name="qcData", tags="*", unitsPerChunk=moreUnits*100000/length(getDataSet(getPlm(this))), moreUnits=1, force=FALSE, verbose=FALSE, ...) {

  # Argument 'path':
  if (is.null(path)) {
    path <- file.path(name, getName(this), "weights", getChipType(getCdf(getDataSet(getPlm(this)))));
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
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  # If residuals already calculated, and if force==FALSE, just return
  # a CelSet with the previous calculations

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filenames <- getNames(getDataSet(this));
  pathname <- sapply(filenames, function(filename){Arguments$getWritablePathname(sprintf("%s,weights.cel", filename), path=path)});
  nbrOfFiles <- length(pathname);
  
  ds <- getDataSet(this);
  ces <- getChipEffects(this);
  paf <- getProbeAffinities(getPlm(this));
  
  nbrOfUnits <- nbrOfUnits(getCdf(ces));

  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

# find number of units to do; first check whether last file in list
# exists

  if (isFile(pathname[nbrOfFiles])) {
    unitsToDo <- findUnitsTodo(QualityAssessmentFile$fromFile(pathname[nbrOfFiles]));
  } else {
    unitsToDo <- 1:nbrOfUnits;
  }

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                    nbrOfChunks, unitsPerChunk);
  head <- 1:unitsPerChunk;
  verbose && enter(verbose, "Extracting unit data");
  count <- 1;
  idxOffset <- as.integer(0);

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

    resFcn <- function(kk) {
      rawData <- rawDataList[[kk]][[1]]$intensities;
      ce <- chipEffectList[[kk]][[1]]$theta[1,];
      pa <- probeAffinityList[[kk]][[1]]$phi;
      fittedData <- outer(pa, ce, FUN="+");
      scale <- median(abs(rawData-fittedData))/0.6745;
      return(matrix(psi.huber((rawData - fittedData)/scale), ncol=nbrOfArrays(ds)));
    }
    
   weightsList <- lapply(head, FUN=resFcn);

# update output files
    
    cdf <- getCellIndices(getCdf(ds), units=units, stratifyBy="pm", ...);
  
    for (kk in seq(pathname)) {
      if (!isFile(pathname[kk])) {
        cdfHeader <- getHeader(getCdf(ds));
        celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=getName(getFile(ds,kk)));
        createCel(pathname[kk], header=celHeader, verbose=less(verbose));
      }
      data <- lapply(weightsList, function(x){nrow <- nrow(x); list(list(intensities=2^x[,kk], stdvs=rep(1, nrow), pixels=rep(1, nrow)))})
      updateCelUnits(pathname[kk], cdf=cdf, data=data);
    }
    
    unitsToDo <- unitsToDo[-head];
    count <- count + 1;
  }

  res <- QualityAssessmentSet$fromFiles(path=path, pattern=",weights.[cC][eE][lL]$");
  setAlias(res, getName(this));
  return(res);
  
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


setMethodS3("plotRle", "QualityAssessmentModel", function(this, subset=NULL, verbose=FALSE, main="RLE", ...) {

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
  verbose && enter(verbose, "Extracting chip effects");
  avg <- getAverageLog(getChipEffects(plm), field="intensities", indices=units, verbose=verbose)
  verbose && exit(verbose);
  medianLE <- getData(avg, indices=units, "intensities")$intensities;
  medianLE <- log2(medianLE);

  verbose && enter(verbose, "Calculating summaries for ", nbrOfArrays(this), " arrays");
  
  # for each file, sweep through and calculate statistics for boxplot.
  boxplotStats <- list();
  for (kk in seq(ces)) {
    chipEffects <- getData(getFile(ces, kk), indices=units, "intensities")$intensities;
    chipEffects <- log2(chipEffects);
    boxplotStats[[kk]] <- boxplot.stats(chipEffects-medianLE);
  }
  rm(avg);
  rm(chipEffects);    # not needed any more

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


##########################################################################
# HISTORY:
# 2007-01-15
# o Renamed to QualityAssessmentModel from QcInfo.
# 2007-01-13
# o Added getWeights().
# 2007-01-06
# o Created.
##########################################################################
