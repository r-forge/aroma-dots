###########################################################################/**
# @RdocClass HetLogAddSnpPlm
#
# @title "The HetLogAddSnpPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "HetLogAddPlm".}
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
setConstructorS3("HetLogAddSnpPlm", function(..., mergeStrands=FALSE, tags="*") {
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


  extend(RmaPlm(..., tags=tags), c("HetLogAddSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


############################################################################
# HISTORY:
# 2007-10-06
# o Created.
############################################################################
