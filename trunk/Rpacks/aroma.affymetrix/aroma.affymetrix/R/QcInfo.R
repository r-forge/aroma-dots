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


setMethodS3("getResiduals", "QcInfo", function(this, unitsPerChunk=moreUnits*100000/length(getDataSet(getPlm(this))), moreUnits=1, ...) {

  # Argument 'unitsPerChunk':
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  ds <- getDataSet(getPlm(this));
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
      str(rawData)
      ce <- chipEffectList[[kk]][[1]]$theta[1,];
      str(ce)
      pa <- probeAffinityList[[kk]][[1]]$phi;
      str(pa)
      fittedData <- outer(pa, ce, FUN="+");
      str(fittedData)
      return(rawData - fittedData);
    }
    
    residualsList <- c(residualsList, lapply(head, FUN=resFcn));
    
    unitsToDo <- unitsToDo[-head];
    count <- count + 1;
  }

  return(residualsList);
  
})

