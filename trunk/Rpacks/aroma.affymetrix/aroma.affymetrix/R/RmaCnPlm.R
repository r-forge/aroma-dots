###########################################################################/**
# @RdocClass RmaCnPlm
#
# @title "The RmaCnPlm class"
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
#   \item{...}{Arguments passed to @see "RmaSnpPlm".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated 
#      seperately, otherwise together.}
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
