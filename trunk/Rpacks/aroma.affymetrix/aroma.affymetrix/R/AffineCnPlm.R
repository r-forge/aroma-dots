###########################################################################/**
# @RdocClass AffineCnPlm
#
# @title "The AffineCnPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffineSnpPlm".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated 
#      seperately, otherwise together.}
#   \item{tags}{A @character @vector of tags.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AffineCnPlm", function(..., combineAlleles=FALSE, tags="*") {
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
      if (combineAlleles)
        tags <- R.utils::insert.default(tags, idx+1, "A+B");
    }
  }


  extend(AffineSnpPlm(..., tags=tags), c("AffineCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})



############################################################################
# HISTORY:
# 2007-01-07
# o Created from MbeiCnPlm.R.
############################################################################
