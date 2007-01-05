###########################################################################/**
# @RdocClass RmaCnPlm
#
# @title "The RmaCnPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaSnpPlm".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated 
#      seperately, otherwise together.}
#   \item{tags}{A @character @vector of tags.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   TO DO.
# }
#
# @author
#
#*/###########################################################################
setConstructorS3("RmaCnPlm", function(..., combineAlleles=FALSE, tags="*") {
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


  extend(RmaSnpPlm(..., tags=tags), c("RmaCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


############################################################################
# HISTORY:
# 2006-09-10
# o Recreated.
############################################################################
