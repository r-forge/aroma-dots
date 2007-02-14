###########################################################################/**
# @RdocClass RmaSnpPlm
#
# @title "The RmaSnpPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaPlm".}
#   \item{mergeStrands}{If @TRUE, the sense and the anti-sense strands are
#      fitted together, otherwise separately.}
#   \item{tags}{A @character @vector of tags.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
#*/###########################################################################
setConstructorS3("RmaSnpPlm", function(..., mergeStrands=FALSE, tags="*") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    idx <- which(tags == "*");
    if (length(idx) > 0) {
      if (!mergeStrands)
        tags <- R.utils::insert.default(tags, idx+1, "+-");
    }
  }


  extend(RmaPlm(..., tags=tags), c("RmaSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})

setMethodS3("getParameterSet", "RmaSnpPlm", function(this, ...) {
  params <- NextMethod("getParameterSet", this, ...);
  params$mergeStrands <- this$mergeStrands;
  params;
}, private=TRUE)



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
