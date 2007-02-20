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
    status <- inferParameters(this, ...)$mergeStrands;

  # Argument 'status':
  status <- Arguments$getLogical(status);

  # Update all chip-effect files
  lapply(this, function(ce) {
    ce$mergeStrands <- status;
  })

  invisible(oldStatus);
})


setMethodS3("inferParameters", "SnpChipEffectSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Infer parameters from stored data in quartet units in CEL set");
  
  # Identify units with quartets
  cdf <- getCdf(this);
  cdfPathname <- getPathname(cdf);
  nbrOfUnits <- readCdfHeader(cdfPathname)$nunits;
  allUnits <- 1:nbrOfUnits;

  ce <- getFile(this, 1);
  cePathname <- getPathname(ce);

  verbose && cat(verbose, "Pathname: ", cePathname);

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
    verbose && cat(verbose, "Scanning units:");
    verbose && str(verbose, units);

    # Infer parameters from 'stdvs'
    stdvs <- readCelUnits(cePathname, units=units, 
                                    readIntensities=FALSE, readStdvs=TRUE);
    # Put quartets by columns
    stdvs <- matrix(unlist(stdvs, use.names=FALSE), nrow=4);
    # Keep only estimated units without NAs
    csums <- colSums(stdvs);
    stdvs <- stdvs[,is.finite(csums) & (csums > 0),drop=FALSE];
    verbose && cat(verbose, "Stdvs quartets:");
    verbose && print(verbose, stdvs[,seq_len(min(ncol(stdvs),6)),drop=FALSE]);
    if (ncol(stdvs) > 0) {
      t <- rowMeans(stdvs);
      if (length(t) > 0) {
        isZero <- isZero(t);
        if (!all(isZero)) {
          mergeStrands <- isZero[3] && isZero[4];
          break;
        }
      }
    }
  } # while(...)

  if (is.na(mergeStrands)) {
    throw("Failed to infer parameter 'mergeStrands' from chip-effect file: ", cePathname);
  }

  res <- list(mergeStrands=mergeStrands);

  verbose && str(verbose, res);
  verbose && exit(verbose);

  res;
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-02-20
# o BUG FIX: inferParameters() would give an error if some estimates were
#   NAs.
# 2007-02-19
# o BUG FIX: If inferParameters() where called on a chip-effect set where
#   some units where not yet estimated, an error would be generated.
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
