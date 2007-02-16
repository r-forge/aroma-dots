###########################################################################/**
# @RdocClass ExonRmaPlm
#
# @title "The ExonRmaPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model part of the Robust Multichip
#  Analysis (RMA) method described in Irizarry et al (2003), as implemented
#  for exon arrays.  The model may be fitted with exons merged into
#  transcripts (all probes fitted together) or on an individual exon basis
#  (probes within an exon treated as a group, but exons fitted separately).
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaPlm".}
#   \item{tags}{A @character @vector of tags.}
#   \item{mergeGroups}{A @logical flag specifying whether to merge exons
#      into transcripts.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#    @see "RmaPlm".
# }
#
# \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
#
# \references{
#  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}. 
#  NAR, 2003, 31, e15.\cr
# }
#*/###########################################################################


setConstructorS3("ExonRmaPlm", function(..., tags="*", mergeGroups=TRUE) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    asteriskTag <- "RMA";
    # Update default tags
    tags[tags == "*"] <- asteriskTag;

    # Split by commas
    tags <- paste(tags, collapse=",");
    tags <- unlist(strsplit(tags, split=","));
  }

  extend(RmaPlm(..., tags=tags), "ExonRmaPlm",
         mergeGroups=mergeGroups
         )
  
})


setMethodS3("getCellIndices", "ExonRmaPlm", function(this, ...) {

# internal functions  
  cdfMergeGroups <- function(groups, ...) {
    nbrOfGroups <- length(groups);

    nbrOfFields <- length(.subset2(groups,1));
    newGroup <- vector("list", nbrOfFields);
    for (ff in seq(length=nbrOfFields)) {
      newGroup[[ff]] <- unlist(lapply(groups, .subset2, ff), use.names=FALSE);
    }
    names(newGroup) <- names(.subset2(groups,1));
    return(list(newGroup));
  }

  cells <- NextMethod("getCellIndices", this, ...);

  # Merge groups?
  if (this$mergeGroups) {
    cells <- applyCdfGroups(cells, cdfMergeGroups);
  }

  cells;
})


setMethodS3("getChipEffectSetClass", "ExonRmaPlm", function(this, ...) {
  ExonChipEffectSet;
}, private=TRUE)


setMethodS3("getChipEffects", "ExonRmaPlm", function(this, ...) {
  ces <- NextMethod("getChipEffects", this, ...);
  setMergeGroups(ces, this$mergeGroups);
  ces;
})

setMethodS3("getProbeAffinities", "ExonRmaPlm", function(this, ..., .class=ExonProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinities", this, ..., .class=.class);
  setMergeGroups(paf, this$mergeGroups);
  paf;
})

setMethodS3("setMergeGroups", "ExonRmaPlm", function(this, ...) {
  ces <- getChipEffects(this);
  setMergeGroups(ces, ...);
  paf <- getProbeAffinities(this);
  setMergeGroups(paf, ...);
})


