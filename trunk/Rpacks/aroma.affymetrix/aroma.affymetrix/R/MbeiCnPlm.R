###########################################################################/**
# @RdocClass MbeiCnPlm
#
# @title "The MbeiCnPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MbeiSnpPlm".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated 
#      seperately, otherwise together.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("MbeiCnPlm", function(..., tags=c("MBEI", ifelse(mergeStrands, "", "+-"), ifelse(combineAlleles, "", "AB")), mergeStrands=FALSE, combineAlleles=FALSE) {
  extend(MbeiSnpPlm(..., tags=tags), c("MbeiCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})



############################################################################
# HISTORY:
# 2006-09-12
# o Recreated.
############################################################################
