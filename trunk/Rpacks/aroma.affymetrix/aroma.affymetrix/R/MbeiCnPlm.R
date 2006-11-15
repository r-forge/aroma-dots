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
#   \item{name}{The name of the PLM, also used as part of the pathname.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("MbeiCnPlm", function(..., combineAlleles=FALSE, name="modelMbeiCnPlm") {
  extend(MbeiSnpPlm(..., name=name), c("MbeiCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})

setMethodS3("getRootPath", "MbeiCnPlm", function(this, ...) {
  "modelMbeiCnPlm";
})


############################################################################
# HISTORY:
# 2006-09-12
# o Recreated.
############################################################################
