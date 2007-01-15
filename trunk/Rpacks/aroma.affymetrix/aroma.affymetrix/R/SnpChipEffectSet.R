###########################################################################/**
# @RdocClass SnpChipEffectSet
#
# @title "The SnpChipEffectSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ChipEffectSet".}
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("SnpChipEffectSet", function(..., mergeStrands=FALSE) {
  this <- extend(ChipEffectSet(...), "SnpChipEffectSet");
  setMergeStrands(this, mergeStrands);
  this;
})

setMethodS3("fromFiles", "SnpChipEffectSet", function(static, ..., mergeStrands="auto") {
  fromFiles.ChipEffectSet(static, ..., mergeStrands=mergeStrands);
}, static=TRUE)



setMethodS3("getAverageFile", "SnpChipEffectSet", function(this, ...) {
  res <- NextMethod(generic="getAverageFile", object=this, ...);
  res$mergeStrands <- getMergeStrands(this);
  res;
})



setMethodS3("getChipEffectFileClass", "SnpChipEffectSet", function(static, ...) {
  SnpChipEffectFile;
}, static=TRUE, private=TRUE)

setMethodS3("getMergeStrands", "SnpChipEffectSet", function(this, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  ce <- getFile(this, 1);
  ce$mergeStrands;
})

setMethodS3("setMergeStrands", "SnpChipEffectSet", function(this, status, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);

  oldStatus <- getMergeStrands(this);

  if (identical(status, "auto"))
    status <- inferParameters(this)$mergeStrands;

  # Argument 'status':
  status <- Arguments$getLogical(status);

  # Update all chip-effect files
  lapply(this, function(ce) {
    ce$mergeStrands <- status;
  })

  invisible(oldStatus);
})


setMethodS3("inferParameters", "SnpChipEffectSet", function(this, ...) {
  # Identify units with quartets
  cdf <- getCdf(this);
  cdfPathname <- getPathname(cdf);
  nbrOfUnits <- readCdfHeader(cdfPathname)$nunits;
  allUnits <- 1:nbrOfUnits;

  ce <- getFile(this, 1);
  cePathname <- getPathname(ce);

  mergeStrands <- NA;
  while(length(allUnits) > 0) {
    uu <- 1:min(1000,length(allUnits));
    units <- allUnits[uu];
    allUnits <- allUnits[-uu];

    # Identify units that are quartets
    unitSizes <- readCdfGroupNames(cdfPathname, units=units);
    names(unitSizes) <- NULL;
    unitSizes <- sapply(unitSizes, FUN=length);
    units <- units[unitSizes == 4];

    # Infer parameters from 'stdvs'
    stdvs <- readCelUnits(cePathname, units=units, 
                                    readIntensities=FALSE, readStdvs=TRUE);
    stdvs <- matrix(unlist(stdvs, use.names=FALSE), nrow=4);
    stdvs <- stdvs[,(colSums(stdvs) > 0),drop=FALSE];
    if (ncol(stdvs) > 0) {
      t <- rowMeans(stdvs);
      isZero <- isZero(t);
      if (!all(isZero)) {
        mergeStrands <- isZero[3] && isZero[4];
        break;
      }
    }
  } # while(...)

  if (is.na(mergeStrands)) {
    throw("Failed to infer parameter 'mergeStrands' from chip-effect file: ", cePathname);
  }

  list(mergeStrands=mergeStrands);
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-01-11
# o Added fromFiles() which now infers 'mergeStrands' from the files.
# o Added inferParameters().
# 2006-11-22
# o Now getAverageFile() finally sets 'mergeStrands'.
# 2006-10-02
# o Added extractSnpQSet() so that we can run crlmm().
# 2006-09-11
# o Created.
############################################################################
