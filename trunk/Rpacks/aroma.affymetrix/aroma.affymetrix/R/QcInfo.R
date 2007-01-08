###########################################################################/**
# @RdocClass QcInfo
#
# @title "The QcInfo class"
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
setConstructorS3("QcInfo", function(plm=NULL, tags="*", ...) {
  
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "QC";
  }

  extend(Object(), "QcInfo",
         .plm = plm,
         .tags = tags );
})

setMethodS3("getPlm", "QcInfo", function(this, ...) {
  this$.plm;
})

setMethodS3("getChipEffects", "QcInfo", function(this, ...) {
  getChipEffects(this$.plm);
})

setMethodS3("getName", "QcInfo", function(this, ...) {
  getName(this$.plm);
})

######################

setMethodS3("as.character", "QcInfo", function(this, ...) {
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

##########################

setMethodS3("nbrOfArrays", "QcInfo", function(this, ...) {
  ces <- getChipEffects(this);
  nbrOfArrays(ces);
})


setMethodS3("getCdf", "QcInfo", function(this, ...) {
  ces <- getChipEffects(this);
  getCdf(ces);
}, protected=TRUE)

setMethodS3("getTags", "QcInfo", function(this, ...) {
  ces <- getChipEffects(this);
  c(getTags(ces), this$.tags);
})

setMethodS3("getFullName", "QcInfo", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})

setMethodS3("getRootPath", "QcInfo", function(this, ...) {
  "qcData";
})

setMethodS3("getDataSet", "QcInfo", function(this, ...) {
  getDataSet(getPlm(this), ...);
})

setMethodS3("getPath", "QcInfo", function(this, ...) {
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
#   Returns an @see "AffymetrixCelSet".
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


setMethodS3("getResiduals", "QcInfo", function(this, path=NULL, tags="*", unitsPerChunk=moreUnits*100000/length(getDataSet(getPlm(this))), moreUnits=1, force=FALSE, ...) {

  # Argument 'path':
  
  
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
  filenames <- basename(getPathnames(getDataSet(this)));
  pathname <- Arguments$getWritablePathname(filename, path=path);
  
  
  ds <- getDataSet(this);
  ces <- getChipEffects(this);
  paf <- getProbeAffinities(getPlm(this));

  nbrOfUnits <- nbrOfUnits(getCdf(ces));

  verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits);

  unitsToDo <- 1:nbrOfUnits;

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n",
                    nbrOfChunks, unitsPerChunk);
  head <- 1:unitsPerChunk;
  verbose && enter(verbose, "Extracting unit data");
  count <- 1;
  idxOffset <- as.integer(0);

  residualsList <- list();
  
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
      return(rawData - fittedData);
    }
    
    residualsList <- c(residualsList, lapply(head, FUN=resFcn));

# update output files
    
    cdf <- getCdfCellIndices(ds, units=units, stratifyBy="pm", ...);
  
    
    
    unitsToDo <- unitsToDo[-head];
    count <- count + 1;
  }

  return(residualsList);
  
})

setMethodS3("plotNuse", "QcInfo", function(this, subset=NULL, verbose=FALSE, ...) {

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
  
  bxp(bxpStats, ...);
  
})


setMethodS3("plotRle", "QcInfo", function(this, subset=NULL, verbose=FALSE, ...) {

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
    boxplotStats[[kk]] <- boxplot.stats(chipEffects/medianLE);
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
  
  bxp(bxpStats, ...);
  
})
