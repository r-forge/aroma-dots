###########################################################################/**
# @RdocClass AffineSnpPlm
#
# @title "The AffineSnpPlm class"
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
#   \item{...}{Arguments passed to @see "AffinePlm".}
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
setConstructorS3("AffineSnpPlm", function(..., mergeStrands=FALSE, tags="*") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Find "replacement" tags
    idx <- which(tags == "*");
    if (length(idx) > 0) {
      if (!mergeStrands)
        tags <- R.utils::insert.default(tags, idx+1, "+-");
    }
  }

  extend(AffinePlm(..., tags=tags), c("AffineSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


############################################################################
# HISTORY:
# 2006-09-11
# o Created from the MbeiSnpPlm.
############################################################################
