###########################################################################/**
# @RdocClass SnpPlm
#
# @title "The SnpPlm interface class"
#
# \description{
#  @classhierarchy
#
#  An @see "Interface" implementing methods special for 
#  @see "ProbeLevelModel"s specific to SNP arrays.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# \section{Requirments}{
#   Classes inheriting from this @see "Interface" must provide the following
#   fields:
#   \itemize{
#    \item{mergeStrands}{A @logical value indicating if strands should be
#       merged or not.}
#   }
# }
#
# @examples "../incl/SnpPlm.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("SnpPlm", function(...) {
  extend(Interface(), "SnpPlm");
})

## setMethodS3("getSubname", "SnpPlm", function(this, ...) {
##   s <- NextMethod("getSubname", this, ...);
##   if (this$mergeStrands) {
##     s <- sprintf("%sStrandless", s);
##   } else {
##     s <- sprintf("%sStrands", s);
##   }
##   s;
## })


setMethodS3("getParameterSet", "SnpPlm", function(this, ...) {
  params <- NextMethod("getParameterSet", this, ...);
  params$mergeStrands <- this$mergeStrands;
  params;
}, private=TRUE)


setMethodS3("getCellIndices", "SnpPlm", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # Merge strands?
  if (this$mergeStrands) {
    cells <- applyCdfGroups(cells, cdfMergeStrands);
  }
  
  cells;
})

setMethodS3("getChipEffectSetClass", "SnpPlm", function(this, ...) {
  SnpChipEffectSet;
}, private=TRUE)


setMethodS3("getChipEffects", "SnpPlm", function(this, ...) {
  ces <- NextMethod("getChipEffects", this, ...);
  setMergeStrands(ces, this$mergeStrands);
  ces;
})

setMethodS3("getProbeAffinities", "SnpPlm", function(this, ..., .class=SnpProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinities", this, ..., .class=.class);
  setMergeStrands(paf, this$mergeStrands);
  paf;
})

setMethodS3("setMergeStrands", "SnpPlm", function(this, ...) {
  ces <- getChipEffects(this);
  setMergeStrands(ces, ...);
  paf <- getProbeAffinities(this);
  setMergeStrands(paf, ...);
})


############################################################################
# HISTORY:
# 2006-09-11
# o The intention is to use SnpPlm as an interface class (that is a class
#   that must not have any fields!) but any class "implementing" this class 
#   (by adding it to its list of classes), will have these methods too.
# o Created from the RmaSnpPlm.
############################################################################
