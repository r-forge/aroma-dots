###########################################################################/**
# @RdocClass RmaSnpPlm
#
# @title "The RmaSnpPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model used in RMA.
#  It can be used to fit the model on a @see "AffymetrixCelSet".
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSet}{An @see "AffymetrixCelSet" object.}
#   \item{path}{The @character string specifying the path to the directory
#      to contain the parameter-estimate files.}
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \section{Model estimates}{
#   The estimated probe affinities are represented by the
#   @see "RmaProbeAffinityFile" class.  
# }
#
# \references{
# }
#*/###########################################################################
setConstructorS3("RmaSnpPlm", function(..., name="modelRmaSnpPlm", mergeStrands=FALSE) {
  this <- extend(RmaPlm(..., name=name), "RmaSnpPlm",
    mergeStrands = mergeStrands
  )
  if (!is.null(getDataSet(this)))
    this <- setup(this);
  this;
})


setMethodS3("setup", "RmaSnpPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  plmCdfMergeStrands <- function(units) {
    if (is.list(units[[1]])) {
      # CDF list structure.
      applyCdfGroups(units, cdfMergeStrands);
    } else {
      # CDF block names: Merge block names for strands
      lapply(units, FUN=function(n) {
        names <- NULL;
        while (n >= 2) {
          name <- paste(unique(n[1:2]), collapse="");
          names <- c(names, name);
          n <- n[-c(1:2)];
        }
        names <- c(names, n);
      });
    }
  } # cdfMergeStrands()
  
  if (this$mergeStrands) {
    # Enforce that the right CDF structure is returned 
    # (for now at least).
    dataSet <- getDataSet(this);
    cdf <- getCdf(dataSet);
    setRestructor(cdf, plmCdfMergeStrands);
    setCdf(dataSet, cdf);
  }

  invisible(this);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-09-11
# o Simple benchmarking [Thinkpad A31]: Fitting 1000 units (with merged 
#   strands) across 22 arrays (100K Xba) takes in total 114 sec, that is,
#   5.1ms/unit/array. 
#   For all 59015 SNPs it takes ~5.0min/array or ~112min/22 arrays.
#   4545 units and 22 arrays: 60s to read all data, 50s to fit the model,
#   30s to store probe affinities, and 120s to store chip-effects. 
#   In total 274s, that is, 2.7ms/unit/array.
#   We are still spending [(60+120)/274 =] 65% on I/O.
#   For all 59015 SNPs it takes ~2.7min/array or ~60min/22 arrays.
# o The fit function now returns all required fields.
# 2006-08-25
# o Created from the corresponding Li & Wong model.
############################################################################
