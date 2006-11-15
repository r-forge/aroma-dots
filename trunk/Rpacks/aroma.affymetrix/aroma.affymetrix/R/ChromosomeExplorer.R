###########################################################################/**
# @RdocClass ChromosomeExplorer
#
# @title "The ChromosomeExplorer class"
#
# \description{
#  @classhierarchy
#
#  This class represents a chromosome explorer.
# }
# 
# @synopsis
#
# \arguments{
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
#*/###########################################################################
setConstructorS3("ChromosomeExplorer", function(ces=NULL, refCe=NULL, ...) {
  if (!is.null(ces)) {
    if (!inherits(ces, "ChipEffectSet")) {
      throw("Argument 'ces' is not a ChipEffectSet object: ", class(ces)[1]);
    }

    if (!inherits(refCe, "ChipEffectFile")) {
      throw("Argument 'refCe' is not a ChipEffectFile object: ", 
                                                             class(refCe)[1]);
    }
  }

  extend(Object(...), "ChromosomeExplorer",
    ces = ces,
    refCe = refCe
  )
})

setMethodS3("getRootPath", "ChromosomeExplorer", function(this, ...) {
  "chromosomeExplorer";
})

setMethodS3("getChipEffects", "ChromosomeExplorer", function(this, ...) {
  this$ces;
})

setMethodS3("getReferenceChipEffect", "ChromosomeExplorer", function(this, ...) {
  this$refCe;
})






############################################################################
# HISTORY:
# 2006-11-02
# o Created.
############################################################################
