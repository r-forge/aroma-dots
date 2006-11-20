###########################################################################/**
# @RdocClass MbeiSnpPlm
#
# @title "The MbeiSnpPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MbeiPlm".}
#   \item{mergeStrands}{If @TRUE, the sense and the anti-sense strands are
#      fitted together, otherwise separately.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("MbeiSnpPlm", function(..., mergeStrands=FALSE, tags="*") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(dataSet)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    idx <- which(tags == "*");
    if (length(idx) > 0) {
      if (!mergeStrands)
        tags <- R.utils::insert.default(tags, idx+1, "+-");
    }
  }


  extend(MbeiPlm(..., tags=tags), c("MbeiSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created from the RmaSnpPlm.
############################################################################
