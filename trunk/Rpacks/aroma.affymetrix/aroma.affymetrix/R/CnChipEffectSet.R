###########################################################################/**
# @RdocClass CnChipEffectSet
#
# @title "The CnChipEffectSet class"
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
#   \item{...}{Arguments passed to @see "SnpChipEffectSet".}
#   \item{combineAlleles}{A @logical indicating if the signals from
#      allele A and allele B are combined or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("CnChipEffectSet", function(..., combineAlleles=FALSE) {
  this <- extend(SnpChipEffectSet(...), "CnChipEffectSet");
  setCombineAlleles(this, combineAlleles);
  this;
})

setMethodS3("fromFiles", "CnChipEffectSet", function(static, ..., combineAlleles="auto") {
  fromFiles.SnpChipEffectSet(static, ..., combineAlleles=combineAlleles);
}, static=TRUE)


setMethodS3("getAverageFile", "CnChipEffectSet", function(this, ...) {
  res <- NextMethod(generic="getAverageFile", object=this, ...);
  res$combineAlleles <- getCombineAlleles(this);
  res;
})

setMethodS3("getChipEffectFileClass", "CnChipEffectSet", function(static, ...) {
  CnChipEffectFile;
}, static=TRUE, private=TRUE)


setMethodS3("getCombineAlleles", "CnChipEffectSet", function(this, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);
  ce <- getFile(this, 1);
  ce$combineAlleles;
})




setMethodS3("setCombineAlleles", "CnChipEffectSet", function(this, status, ...) {
  if (nbrOfFiles(this) == 0)
    return(FALSE);

  ce <- getFile(this, 1);
  oldStatus <- ce$ombineAlleles;

  if (identical(status, "auto"))
    status <- inferParameters(this)$combineAlleles;

  status <- Arguments$getLogical(status);
  lapply(this, FUN=function(ce) {
    ce$combineAlleles <- status;
  })
  invisible(oldStatus);
})


setMethodS3("inferParameters", "CnChipEffectSet", function(this, ...) {
  # Identify units with quartets
  cdf <- getCdf(this);
  cdfPathname <- getPathname(cdf);
  nbrOfUnits <- readCdfHeader(cdfPathname)$nunits;
  allUnits <- 1:nbrOfUnits;

  ce <- getFile(this, 1);
  cePathname <- getPathname(ce);

  mergeStrands <- combineAlleles <- NA;
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
        combineAlleles <- isZero[2] && isZero[4];
        mergeStrands <- isZero[3] && isZero[4];
        break;
      }
    }
  } # while(...)

  if (is.na(mergeStrands) || is.na(combineAlleles)) {
    throw("Failed to infer parameters 'mergeStrands' and 'combineAlleles' from chip-effect file: ", cePathname);
  }

  list(combineAlleles=combineAlleles, mergeStrands=mergeStrands);
}, private=TRUE)




############################################################################
# HISTORY:
# 2007-01-11
# o Added fromFiles() which now infers 'combineAlleles' from the files.
# o Added inferParameters().
# 2006-11-22
# o Now getAverageFile() finally sets 'combineAlleles'.
# 2006-09-11
# o Created.
############################################################################
